import os
import sys
import threading
import json
import uuid
import sqlite3

script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_directory, '..'))  # Add the parent directory to the path

from werkzeug.security import generate_password_hash, check_password_hash
from flask_sqlalchemy import SQLAlchemy
from scripts.upload.upload_script import send_sparql_query
from scripts.upload.upload_script import upload_binary_to_existing_docker
from scripts.mapping.mixmapping import mappingmixture
from scripts.mapping.mappingscript import placeholderreplacement
from scripts.rawdataextraction.mixdesign_metadata_extraction import mix_metadata
from datetime import timedelta, datetime
from flask import Flask, request, render_template, redirect, url_for, session, jsonify



app = Flask(__name__)

# Load config.json
with open('config.json', 'r') as file:
    config = json.load(file)
    app.config['SECRET_KEY'] = config['SECRET_KEY']

# clear users.db before deployment !!!
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
app.config['SESSION_PERMANENT'] = True
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=3650)  # 10 years
app.config['SESSION_TYPE'] = "filesystem"

db = SQLAlchemy(app)

# Upload Database
upload_db = 'upload.db'


# Create password db
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    password_hash = db.Column(db.String(128), nullable=False)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)


with app.app_context():
    db.create_all()


# main page
@app.route('/')
def index():
    # go to welcome page (no login required)
    return render_template('welcome.html')


# login page
@app.route('/login', methods=['GET', 'POST'])
def login():
    if 'username' in session:
        return redirect(url_for('query'))
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        user = User.query.filter_by(username=username).first()
        if user and user.check_password(password):  # check password
            session['username'] = user.username
            # set session for user
            return redirect(url_for('query_page'))
        else:
            error_message = 'Incorrect Username or Password'
            return render_template('login.html', errormessage=error_message)
    return render_template('login.html')


# Sign up page
@app.route('/signup', methods=['GET', 'POST'])
def signup():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        # Check if the username already exists
        existing_user = User.query.filter_by(username=username).first()
        if existing_user:
            # Username already taken, show an error message
            error_message = "Username already exists. Please choose a different one."
            return render_template('signup.html', errormessage=error_message)

        # Create a new user
        new_user = User(username=username)
        new_user.set_password(password)

        # Add to the database
        db.session.add(new_user)
        db.session.commit()
        return redirect(url_for('login'))

    return render_template('signup.html')


# Send Query to Ontodocker
@app.route('/queryexec', methods=['POST'])
def execute_sparql_query():
    if 'username' in session:
        # Hier extrahieren wir die Daten aus dem POST-Request
        results = send_sparql_query(request.form['query'], config)
        return jsonify({
            'message': results
        })


# Query page
@app.route('/query')
def query_page():
    if 'username' in session:
        return render_template('queryPage.html')
    else:
        return redirect(url_for('login'))


# Query page (simple version)
@app.route('/query-simple')
def query_page_simple():
    if 'username' in session:
        return render_template('queryPage-simple.html')
    else:
        return redirect(url_for('login'))


# Upload page
@app.route('/upload')
def upload_page():
    if 'username' in session:
        return render_template('uploadForm.html')
    else:
        return redirect(url_for('login'))


# Logout
@app.route('/logout')
def logout():
    # delete username out of session
    session.pop('username', None)

    return redirect(url_for('login'))


# query mixture (in upload page)
@app.route('/search-mixture', methods=['POST'])
def search_mixture():
    if 'username' in session:
        data = request.json
        mixture_name = data['mixtureName']
        # Search, if the mixture exists in the database
        query = f'''
        SELECT ?s WHERE {{
            ?s <https://w3id.org/pmd/co/value> "{mixture_name}" .
            ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/pmd/co/ProvidedIdentifier> .
          }}
        '''
        value_of_i = 0
        results = send_sparql_query(query, config)
        if results and results.get("results", {}).get("bindings", []):
            # If found, find the ID
            for result in results.get('results', {}).get('bindings', []):
                # Extrahieren des Wertes für 's'
                value = result.get('s', {}).get('value', '')
                # Überprüfen, ob "humanreadableID" im Wert enthalten ist
                if "humanreadableID" in value:
                    query = f'''
                        SELECT ?i WHERE {{
                            ?s <https://w3id.org/pmd/co/value> "{mixture_name}" . ?s <http://www.w3.org/1999/02/22
                            -rdf-syntax-ns#type> <https://w3id.org/pmd/co/ProvidedIdentifier> . ?a 
                            <http://purl.org/spar/datacite/hasIdentifier> ?s . ?a 
                            <http://purl.org/spar/datacite/hasIdentifier> ?b . FILTER(?s != ?b) ?b 
                            <https://w3id.org/pmd/co/value> ?i }} '''
                    results = send_sparql_query(query, config)
                    if results and results.get("results", {}).get("bindings", []):
                        # Extrahieren des Wertes von `?i`
                        for binding in results['results']['bindings']:
                            value_of_i = binding['i']['value']
                        if value_of_i:
                            return jsonify({'message': f'Mischung {mixture_name} erfolgreich gefunden!',
                                            'mixtureID': value_of_i})

                else:
                    value_of_i = mixture_name
            return jsonify({'message': f'Mischung {mixture_name} erfolgreich gefunden!', 'mixtureID': value_of_i})
        else:
            # Keine Ergebnisse, Mischung nicht gefunden
            return jsonify({'message': 'Mischung nicht gefunden.'})


def get_db_connection():
    try:
        conn = sqlite3.connect(upload_db)
        conn.row_factory = sqlite3.Row
    except sqlite3.Error as e:
        print(f"Database connection error: {e}")
        conn = None
    return conn


def init_db():
    with sqlite3.connect(upload_db) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS uploads (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                Unique_ID TEXT,
                Mixture_ID TEXT,
                user TEXT,
                filetype TEXT,
                filename TEXT,
                type TEXT,
                blob BLOB,
                Json BLOB,
                ttl BLOB,
                Json_Specimen BLOB,
                ttl_Specimen BLOB,
                UploadDate TEXT,
                Mapped INTEGER,
                Error INTEGER
            );
        ''')
        conn.commit()


# Mapping function
def async_function(unique_id):

    # Lookup table for path
    paths = {'EModule': 'cpto/EModuleOntology_KG_Template.ttl',
             'Specimen': 'cpto/Specimen_KG_Template.ttl'}

    def add_data(rowname, data):
        # Verbindung zur Datenbank herstellen
        conn = get_db_connection()
        cursor = conn.cursor()
        # UPDATE-Anweisung vorbereiten, um den json-Wert in der Zeile mit der gegebenen uniqueID zu aktualisieren
        query = f"UPDATE uploads SET {rowname} = ? WHERE Unique_ID = ?"
        try:
            cursor.execute(query, (data, unique_id))
            conn.commit()
            print(f"{rowname}-Wert für {unique_id} erfolgreich aktualisiert.")
            success = 1
        except Exception as ex:
            print(f"Fehler beim Aktualisieren des {rowname}-Werts für {unique_id}: {ex}")
            conn.rollback()
            success = 0
        finally:
            # Verbindung schließen
            conn.close()

        return success

    def get_data():
        # Verbindung zur Datenbank herstellen und Daten auslesen
        conn = get_db_connection()
        cursor = conn.cursor()

        # SQL-Abfrage vorbereiten, um alle Daten aus der Zeile mit der gegebenen uniqueID zu erhalten
        query = "SELECT * FROM uploads WHERE Unique_ID = ?"
        cursor.execute(query, (unique_id,))

        # Ergebnis der Abfrage abrufen
        rowdata = cursor.fetchone()
        # Verbindung schließen
        conn.close()

        return rowdata

    def upload_to_docker(data):
        success = upload_binary_to_existing_docker(data, config)
        if success != 0:
            add_data('Error', 1)

    # Imitiere eine lang laufende Aufgabe
    print(f"Starte die Verarbeitung für {unique_id}")

    # Ergebnis der Abfrage abrufen
    row = get_data()

    # Test ob bereits gemapped
    if row['Mapped'] != 0:
        print(f'Error already mapped with {unique_id}')
        return

    # wenn eine json datei hochgeladen wurde
    if row['filetype'] == 'json':
        json_data = json.loads(row['blob'].decode('utf-8'))
        json_data['ID'] = unique_id
        if row['type'] != 'Mixture':
            json_data['mixtureID'] = row['Mixture_ID']
        add_data('Json', json.dumps(json_data).encode('utf-8'))
    else:
        if row['type'] == 'Mixture':
            json_data = mix_metadata(row['blob'], row['filename'])
            json_data['ID'] = row['unique_id']
            print(json_data)
            add_data('Json', json.dumps(json_data).encode('utf-8'))
        # Raw data extraction
        print("RawData")

    # Set ID and mixtureID in json!
    # Mapping
    # Aktualisierten Daten abfragen
    row = get_data()

    if row['Json'] == '':
        print(f"No Json for {unique_id}")
        add_data('Error', 1)

    if row['type'] == 'Mixture':
        try:
            ttl_blob = mappingmixture(row['Json'])
            status = add_data('ttl', ttl_blob)
            # if error set error
            if status != 1:
                add_data('Error', 1)
                return
            else:
                upload_to_docker(ttl_blob)
        except Exception as e:
            print(f'Error in Mixture mapping: {e}')
    else:
        try:
            ttl_blob = placeholderreplacement(paths[row['type']], row['Json'])
            status = add_data('ttl', ttl_blob)
            if status != 1:
                add_data('Error', 1)
                return
            else:
                upload_to_docker(ttl_blob)
        except Exception as e:
            print(f'Error in Placeholderreplacement: {e}')
            add_data('Error', 1)
            return
        if row['Json_Specimen'] != '':
            try:
                ttl_blob = placeholderreplacement(paths['Specimen'], row['Json_Specimen'])
                status = add_data('ttl_Specimen', ttl_blob)
                if status != 1:
                    add_data('Error', 1)
                    return
                else:
                    upload_to_docker(ttl_blob)
            except Exception as e:
                print(f'Error in Placeholderreplacement: {e}')
                add_data('Error', 1)
                return

    print(f"Verarbeitung für {unique_id} abgeschlossen.")
    add_data('Mapped', 1)


# handle uploaded data
@app.route('/dataUpload', methods=['POST'])
def data_upload():
    init_db()
    file_types = ['xlsx', 'xls', 'csv', 'dat', 'txt', 'json']
    if 'username' not in session:
        return jsonify({'error': 'Nicht angemeldet'}), 403

    if 'file' not in request.files:
        return jsonify({'error': 'Keine Datei gefunden'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'Keine Datei ausgewählt'}), 400

    if 'type' not in request.form:
        return jsonify({'error': 'Kein Typ angegeben'}), 400

    if 'Mixture_ID' not in request.form:
        return jsonify({'error': 'Keine Mixture ID angegeben'}), 400

    type = request.form['type']
    mixtureID = str(request.form['Mixture_ID'])
    user = session['username']  # Benutzernamen aus der Session holen

    file_name = file.filename
    # Extrahieren der Dateiendung
    _, file_extension = os.path.splitext(file.filename)

    filetype = file_extension.lstrip('.')

    if filetype not in file_types:
        return jsonify({'error': 'Unsupported Type'}), 400

    # Datei als BLOB speichern
    file_blob = file.read()

    # Aktuelle Zeit
    uploaddate = datetime.now().isoformat()

    if type == 'Mixture':
        unique_id = mixtureID
    else:
        unique_id = str(uuid.uuid4())

    # Verbindung zur Datenbank herstellen und die Daten speichern
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute('INSERT INTO uploads (user, filetype, filename, type, blob, Mixture_ID, Unique_ID, UploadDate, '
                   'Mapped, Error) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                   (user, filetype, file_name, type, file_blob, mixtureID, unique_id, uploaddate, 0, 0))
    conn.commit()
    conn.close()

    # Starten der asynchronen Funktion
    thread = threading.Thread(target=async_function, args=(unique_id,))
    thread.start()

    return jsonify({'message': f'Datei {file.filename} und Typ {type} erfolgreich hochgeladen und gespeichert!',
                    'uniqueID': unique_id}), 200


if __name__ == '__main__':
    app.run(debug=True)
