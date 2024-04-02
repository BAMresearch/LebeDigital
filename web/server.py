import os
import sys
# Get the directory of the current script
script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_directory, '..'))  # Add the parent directory to the path

import uuid
import sqlite3
from datetime import timedelta, datetime
from flask import Flask, request, render_template, redirect, url_for, flash, session, jsonify
from flask_bootstrap import Bootstrap
from werkzeug.security import generate_password_hash, check_password_hash
from flask_sqlalchemy import SQLAlchemy
from scripts.upload.upload_script import send_sparql_query

app = Flask(__name__)

# clear users.db before deployment !!!
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
# change before deployment !!!!!
app.config['SECRET_KEY'] = '9dddba0645ea2c32c78bd5c3409ec5b5c55efc9f9aeb16b750dfd72237e77c8fe52ca1655048fdbff5f5707716670a0cf3e0a01f59a704a1194a59ddc41df693'
# !!!!
app.config['SESSION_PERMANENT'] = True
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=3650)  # 10 years
app.config['SESSION_TYPE'] = "filesystem"

db = SQLAlchemy(app)

# Upload Database
upload_db = 'upload.db'

main_path = ''


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


@app.route('/login', methods=['GET', 'POST'])
def login():
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


#Sign up page
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
        return render_template('login.html')

    return render_template('signup.html')



# Query page
@app.route('/query')
def query_page():
    if 'username' in session:
        return render_template('queryPage.html')
    else:
        return redirect(url_for('login'))


# Upload page
@app.route('/upload')
def upload_page():
    if 'username' in session:
        return render_template('uploadForm.html')
    else:
        return redirect(url_for('login'))


@app.route('/logout')
def logout():
    # delete username out of session
    session.pop('username', None)

    return redirect(url_for('login'))


# main page
@app.route('/')
def index():
    # standard is query (no login required)
    return redirect(url_for('query_page'))


# query mixture
@app.route('/search', methods=['POST'])
def search():
    data = request.json
    mixture_name = data['mixtureName']
    # Search, if the mixture exists in the database
    query = f'''
    SELECT ?s WHERE {{
        ?s <https://w3id.org/pmd/co/value> "{mixture_name}" .
        ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/pmd/co/ProvidedIdentifier> .
      }}
    '''
    print(query)
    value_of_i = 0
    results = send_sparql_query(query)
    if results and results.get("results", {}).get("bindings", []):
        # If found, find the ID
        for result in results.get('results', {}).get('bindings', []):
            # Extrahieren des Wertes für 's'
            value = result.get('s', {}).get('value', '')
            # Überprüfen, ob "humanreadableID" im Wert enthalten ist
            if "humanreadableID" in value:
                query = f'''
                    SELECT ?i WHERE {{
                        ?s <https://w3id.org/pmd/co/value> "{mixture_name}" .
                        ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/pmd/co/ProvidedIdentifier> .
                        ?a <http://purl.org/spar/datacite/hasIdentifier> ?s .
                        ?a <http://purl.org/spar/datacite/hasIdentifier> ?b .
                        FILTER(?s != ?b)
                        ?b <https://w3id.org/pmd/co/value> ?i
                  }}
                '''
                results = send_sparql_query(query)
                if results and results.get("results", {}).get("bindings", []):
                    # Extrahieren des Wertes von `?i`
                    for binding in results['results']['bindings']:
                        value_of_i = binding['i']['value']
                    if value_of_i:
                        return jsonify({'message': f'Mischung {mixture_name} erfolgreich gefunden!', 'mixtureID': value_of_i})

            else:
                value_of_i = mixture_name
        return jsonify({'message': f'Mischung {mixture_name} erfolgreich gefunden!', 'mixtureID': value_of_i})
    else:
        # Keine Ergebnisse, Mischung nicht gefunden
        return jsonify({'message': 'Mischung nicht gefunden.'})


def get_db_connection():
    conn = sqlite3.connect(upload_db)
    conn.row_factory = sqlite3.Row  # Ermöglicht den Zugriff auf die Daten per Index und per Name
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


# handle uploaded data
@app.route('/dataUpload', methods=['POST'])
def data_upload():
    init_db()
    file_types = ['xlsx', 'xls', 'csv', 'dat', 'txt']
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
    cursor.execute('INSERT INTO uploads (user, filetype, type, blob, Mixture_ID, Unique_ID, UploadDate, Mapped, Error) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
                   (user, filetype, type, file_blob, mixtureID, unique_id, uploaddate, 0, 0))
    conn.commit()
    conn.close()

    return jsonify({'message': f'Datei {file.filename} und Typ {type} erfolgreich hochgeladen und gespeichert!'}), 200


if __name__ == '__main__':
    app.run(debug=True)
