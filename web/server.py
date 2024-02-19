import os
import sqlite3
from datetime import timedelta
from flask import Flask, request, render_template, redirect, url_for, flash, session, jsonify
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
            return redirect(url_for('upload_page'))
        else:
            flash('Falscher Benutzername oder Passwort')
    return render_template('login.html')


# Query page
@app.route('/query')
def query_page():
    return render_template('query_page.html')


# Query page
@app.route('/upload')
def upload_page():
    if 'username' in session:
        return render_template('upload_form.html')
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
                user TEXT NOT NULL,
                filetype TEXT NOT NULL,
                type TEXT NOT NULL,
                blob BLOB NOT NULL
            );
        ''')
        conn.commit()


# handle uploaded data
@app.route('/dataUpload', methods=['POST'])
def data_upload():
    init_db()
    file_types = ['xlsx', 'xls', 'csv', 'dat', 'txt']
    # Session-Beispiel (stellen Sie sicher, dass Sie den Benutzernamen in der Session setzen)
    if 'username' not in session:
        return jsonify({'error': 'Nicht angemeldet'}), 403

    if 'file' not in request.files:
        return jsonify({'error': 'Keine Datei gefunden'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'Keine Datei ausgewählt'}), 400

    if 'type' not in request.form:
        return jsonify({'error': 'Kein Typ angegeben'}), 400
    type = request.form['type']
    user = session['username']  # Benutzernamen aus der Session holen
    # Extrahieren der Dateiendung
    _, file_extension = os.path.splitext(file.filename)

    filetype = file_extension.lstrip('.')

    if filetype not in file_types:
        return jsonify({'error': 'Unsupported Type'}), 400

    # Datei als BLOB speichern
    file_blob = file.read()

    # Verbindung zur Datenbank herstellen und die Daten speichern
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute('INSERT INTO uploads (user, filetype, type, blob) VALUES (?, ?, ?, ?)',
                   (user, filetype, type, file_blob))
    conn.commit()
    conn.close()

    return jsonify({'message': f'Datei {file.filename} und Typ {type} erfolgreich hochgeladen und gespeichert!'}), 200


if __name__ == '__main__':
    app.run(debug=True)
