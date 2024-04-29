import os
import sys

script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_directory, '..'))  # Add the parent directory to the path

import threading
import json
import uuid
import sqlite3
from werkzeug.security import generate_password_hash, check_password_hash
from flask_sqlalchemy import SQLAlchemy
from scripts.upload.upload_script import send_sparql_query, upload_binary_to_existing_docker, clear_dataset
from scripts.mapping.mixmapping import mappingmixture
from scripts.mapping.mappingscript import placeholderreplacement
from scripts.rawdataextraction.emodul_xml_to_json import xml_to_json
from scripts.rawdataextraction.mixdesign_metadata_extraction import mix_metadata
from datetime import timedelta, datetime
from flask import Flask, request, render_template, redirect, url_for, session, jsonify, abort, send_file
from loguru import logger
from pathlib import Path
import requests
from werkzeug.datastructures import FileStorage
from io import BytesIO
import zipfile


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


# Pfad zur Basis des Projekts bestimmen
baseDir = Path(__file__).parent
logDir = os.path.join(baseDir, "logs")

# Stelle sicher, dass das Verzeichnis für die Logs existiert
os.makedirs(logDir, exist_ok=True)

# Pfade für die Log-Dateien
debugLogPath = os.path.join(logDir, "debug_{time}.log")
infoLogPath = os.path.join(logDir, "info_{time}.log")

# Konfiguriere die globalen Logger-Einstellungen
logger.configure(handlers=[
    {"sink": debugLogPath, "level": "DEBUG", "rotation": "10 MB", "retention": "10 days"},
    {"sink": infoLogPath, "level": "INFO", "rotation": "10 MB"}
])


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
    #if 'username' in session:
    #    return redirect(url_for('query_page'))
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


# Upload page
@app.route('/admin')
def admin_page():
    if session['username'] == 'admin':
        return render_template('admin.html')
    else:
        return redirect(url_for('login'))


# Logout
@app.route('/logout')
def logout():
    # delete username out of session
    session.pop('username', None)

    return redirect(url_for('login'))


@app.route('/adminData', methods=['POST', 'GET'])
def get_admin_data():
    if session['username'] != 'admin':
        return jsonify({"error": "Unauthorized"}), 401

    if request.method == 'GET':
        # Logik für GET-Anfragen
        users = User.query.all()
        usernames = [user.username for user in users]
        return jsonify({"config": config, "users": usernames})
    elif request.method == 'POST':
        # Logik für POST-Anfragen
        data = request.json
        print(data)
        if data.get("clearData"):
            if data["clearData"]:
                clear_dataset(config)
                print("Ontodocker cleared.")
                # Verbindung zur Datenbank herstellen
                conn = get_db_connection()
                cursor = conn.cursor()
                query = f"DELETE FROM uploads"
                try:
                    cursor.execute(query)
                    conn.commit()
                    print("Datenbank gelöscht.")
                except Exception as ex:
                    print(f"Fehler beim Löschen der Tabelle uploads: {ex}")
                    conn.rollback()
                finally:
                    # Verbindung schließen
                    conn.close()
        else:
            config['ontodocker_url'] = data["config"]["ontodocker_url"]
            config["DOCKER_TOKEN"] = data["config"]["DOCKER_TOKEN"]
            config["triplestore_server"] = data["config"]["triplestore_server"]
            config["dataset_name"] = data["config"]["dataset_name"]

            with open('config.json', 'w') as file:
                json.dump(config, file, indent=4)

            for entry in data["users"]:
                # Suche den Benutzer anhand des Benutzernamens
                user = User.query.filter_by(username=entry).first()

                if user:
                    # Wenn der Benutzer gefunden wurde, lösche ihn
                    db.session.delete(user)
                    db.session.commit()
                else:
                    # Wenn kein Benutzer gefunden wurde, sende eine Fehlermeldung
                    return jsonify({'error': 'User not found'}), 404

        return jsonify({"message": "Data received"})
    else:
        return jsonify({"error": "Method Not Allowed"}), 405


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
        cursor.execute('''
                    CREATE TABLE IF NOT EXISTS uidlookup (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        Unique_ID TEXT,
                        Name TEXT
                    );
                ''')
        conn.commit()


# Mapping function
def async_function(unique_id):

    # Lookup table for path
    paths = {'EModule': '../cpto/EModuleOntology_KG_Template.ttl',
             'Specimen': '../cpto/Specimen_KG_Template.ttl'}

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

    def get_data(uid=None):
        # Verbindung zur Datenbank herstellen und Daten auslesen
        conn = get_db_connection()
        cursor = conn.cursor()

        if uid is None:
            uid = unique_id

        # SQL-Abfrage vorbereiten, um alle Daten aus der Zeile mit der gegebenen uniqueID zu erhalten
        query = "SELECT * FROM uploads WHERE Unique_ID = ?"
        cursor.execute(query, (uid,))

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
        elif row['type'] == 'EModule':
            if row['filetype'] == 'xml':
                mix_data = get_data(row['Mixture_ID'])
                print(mix_data)
                json_data = xml_to_json(row['blob'], mix_data['Json'])
                print(json_data)
                emodule_json = json_data[0]
                specimen_json = json_data[1]
                emodule_json['ID'] = row['unique_id']
                specimen_json['ID'] = row['unique_id']
                emodule_json['SpecimenID'] = row['unique_id']
                print(emodule_json)
                print(specimen_json)
                add_data('Json', json.dumps(emodule_json).encode('utf-8'))
                add_data('Json_Specimen', json.dumps(specimen_json).encode('utf-8'))
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
    file_types = ['xlsx', 'xls', 'csv', 'dat', 'txt', 'json', 'xml']
    if 'username' not in session:
        return jsonify({'error': 'Nicht angemeldet'}), 403

    print(request.files)
    if 'file' in request.files and request.files['file'].filename != '':
        file = request.files['file']
        file_name = file.filename
    elif 'url' in request.form and request.form['url'] != '':
        url = request.form['url']
        response = requests.get(url)
        file = FileStorage(BytesIO(response.content), filename=url.split('/')[-1])
        file_name = url
    else:
        return jsonify({'error': 'Keine Datei oder URL gefunden'}), 400

    if 'type' not in request.form:
        return jsonify({'error': 'Kein Typ angegeben'}), 400

    if 'Mixture_ID' not in request.form:
        return jsonify({'error': 'Keine Mixture ID angegeben'}), 400

    type = request.form['type']
    mixtureID = str(request.form['Mixture_ID'])
    user = session['username']  # Benutzernamen aus der Session holen

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
    cursor.execute('INSERT INTO uidlookup (Unique_ID, Name) VALUES (?, ?)',
                   (unique_id, file_name.split(".")[0]))
    conn.commit()
    conn.close()

    # Starten der asynchronen Funktion
    thread = threading.Thread(target=async_function, args=(unique_id,))
    thread.start()

    return jsonify({'message': f'Datei {file.filename} und Typ {type} erfolgreich hochgeladen und gespeichert!',
                    'uniqueID': unique_id}), 200


@app.route('/rawdownload')
def raw_download():

    temp_directory = "temp/"

    if not os.path.exists(temp_directory):
        os.makedirs(temp_directory)

    if not os.path.exists("zip/"):
        os.makedirs("zip/")

    def zip_directory(folder_path, output_zip_path):
        # Erstelle eine ZIP-Datei und füge Dateien hinzu
        with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Durchlaufe alle Ordner und Dateien im Verzeichnis
            for foldername, subfolders, filenames in os.walk(folder_path):
                for filename in filenames:
                    # Erstelle den absoluten Pfad zur Datei
                    file_path = os.path.join(foldername, filename)
                    # Bestimme den relativen Pfad zur Datei im Verzeichnis für die ZIP-Datei
                    arcname = os.path.relpath(file_path, start=folder_path)
                    # Füge die Datei zur ZIP-Datei hinzu
                    zipf.write(file_path, arcname=arcname)

        logger.info(f"ZIP-Datei wurde erstellt: {output_zip_path}")

    def get_data(uid, row="Unique_ID", table="uploads"):
        # Verbindung zur Datenbank herstellen und Daten auslesen
        conn = get_db_connection()
        cursor = conn.cursor()

        # SQL-Abfrage vorbereiten, um alle Daten aus der Zeile mit der gegebenen uniqueID zu erhalten
        query = f"SELECT * FROM {table} WHERE {row} = ?"
        cursor.execute(query, (uid,))

        # Ergebnis der Abfrage abrufen
        rowdata = cursor.fetchone()
        # Verbindung schließen
        conn.close()

        return rowdata

    # Hole die ID aus den Query-Parametern
    mixture_id = request.args.get('id')

    if not mixture_id:
        abort(400, description="No file ID provided.")

    try:
        uuid_obj = uuid.UUID(mixture_id, version=4)
        if not uuid_obj.version == 4:
            abort(400, description="Wrong ID format")
    except ValueError:
        row = get_data(mixture_id, "Name", "uidlookup")
        mixture_id = row["Unique_ID"]
    if mixture_id:
        row = get_data(mixture_id)

        if row["blob"]:
            # Pfad, unter dem die extrahierten Dateien gespeichert werden
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}.{row['filetype']}")

            # Schreibe BLOB-Daten in eine Datei
            with open(output_path, 'wb') as file:
                file.write(row['blob'])

        if row["Json"]:
            # Pfad, unter dem die extrahierten Dateien gespeichert werden
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}.json")

            # Schreibe BLOB-Daten in eine Datei
            with open(output_path, 'wb') as file:
                file.write(row['Json'])

        if row["ttl"]:
            # Pfad, unter dem die extrahierten Dateien gespeichert werden
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}.ttl")

            # Schreibe BLOB-Daten in eine Datei
            with open(output_path, 'wb') as file:
                file.write(row['ttl'])

        if row["Json_Specimen"]:
            # Pfad, unter dem die extrahierten Dateien gespeichert werden
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}_Specimen.json")

            # Schreibe BLOB-Daten in eine Datei
            with open(output_path, 'wb') as file:
                file.write(row["Json_Specimen"])

        if row["ttl_Specimen"]:
            # Pfad, unter dem die extrahierten Dateien gespeichert werden
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}_Specimen.ttl")

            # Schreibe BLOB-Daten in eine Datei
            with open(output_path, 'wb') as file:
                file.write(row["ttl_Specimen"])

        if row["Mixture_ID"] and row["Mixture_ID"] != row["Unique_ID"]:
            row2 = get_data(row["Mixture_ID"])

            # Pfad, unter dem die extrahierten Dateien gespeichert werden
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}_Mixture.{row2['filetype']}")

            # Schreibe BLOB-Daten in eine Datei
            with open(output_path, 'wb') as file:
                file.write(row2["blob"])

        zip_output_path = f"zip/{row['filename'].split('.')[0]}.zip"
        zip_directory(temp_directory, zip_output_path)

        try:
            # Sende die Datei zum Download, stellt sicher, dass as_attachment=True gesetzt ist
            return send_file(zip_output_path, as_attachment=True)
        except Exception as e:
            logger.warning(f"Error in Zip: {e}")
            # Wenn etwas schiefgeht, z.B. Datei nicht gefunden
            abort(404, description="File not found.")
        #finally:
        #    os.remove(zip_output_path)
        #    os.remove(temp_directory)
    else:
        # Wenn keine ID angegeben ist, sende einen 400 Bad Request Fehler
        abort(400, description="No file ID provided.")

if __name__ == '__main__':
    app.run(debug=True)
