import os
import sys

script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_directory, '..'))  # Add the parent directory to the path

import threading
import json
import uuid
import sqlite3
import time
from werkzeug.security import generate_password_hash, check_password_hash
from flask_sqlalchemy import SQLAlchemy
from scripts.upload.upload_script import send_sparql_query, upload_binary_to_existing_docker, clear_dataset, delete_specific_triples_from_endpoint
from scripts.mapping.mixmapping import mappingmixture
from scripts.mapping.mappingscript import placeholderreplacement
from scripts.rawdataextraction.emodul_xml_to_json import xml_to_json
from scripts.rawdataextraction.mixdesign_metadata_extraction import mix_metadata
from scripts.rawdataextraction.ComSt_generate_processed_data import processed_rawdata
#from scripts.rawdataextraction.ComSt_metadata_extraction import extract_metadata_ComSt
from datetime import timedelta, datetime
from flask import Flask, request, render_template, redirect, url_for, session, jsonify, abort, send_file
from loguru import logger
from pathlib import Path
import requests
from werkzeug.datastructures import FileStorage
from io import BytesIO
import zipfile
import shutil


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


# path to the logs directory
baseDir = Path(__file__).parent
logDir = os.path.join(baseDir, "logs")

# make sure the logs directory exists
os.makedirs(logDir, exist_ok=True)

# path to the log files
debugLogPath = os.path.join(logDir, "debug_{time}.log")
infoLogPath = os.path.join(logDir, "info_{time}.log")

# configure the logger
logger.configure(handlers=[
    {"sink": debugLogPath, "level": "DEBUG", "rotation": "10 MB", "retention": "5 days"},
    {"sink": infoLogPath, "level": "INFO", "rotation": "10 MB"}
])

# delete empty log files older than two days
def delete_empty_logs(directory):
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # check if the file is empty
        if os.path.isfile(filepath) and os.path.getsize(filepath) == 0:
            # check if the file is older than two days
            if time.time() - os.path.getmtime(filepath) > 2 * 24 * 60 * 60:
                os.remove(filepath)

delete_empty_logs(logDir)

# Create password db
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    password_hash = db.Column(db.String(128), nullable=False)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)


# extracts all json information from the database in two lists
def create_json_file():
    '''
    This function extracts all json information from the database.
    The first list contains the json information of the uploads table.
    The second list contains the json_specimen information of the uploads table.
    
    :return: One Json Object with keys for the types (Mixture eg) and two lists as values containing lists of all info of json and json_specimen.
    '''

    typeliste = ['Mixture', 'EModule']
    json_full = {}

    for entry in typeliste:
        conn = get_db_connection()
        cursor = conn.cursor()
        if entry == 'Mixture':
            query = "SELECT Json FROM uploads WHERE type = ?"
            cursor.execute(query, (entry,))
            # Fetch the results of the query
            rows = cursor.fetchall()
            conn.close()
            # Convert the BLOB data to JSON dicts
            json_full[entry] = [json.loads(row['Json'].decode('utf-8')) for row in rows]

        else:
            query = "SELECT u1.Json, u1.Json_specimen, u2.Json FROM uploads u1 INNER JOIN uploads u2 ON u1.Mixture_ID = u2.Unique_ID WHERE u1.type = ?"
            cursor.execute(query, (entry,))
            # Fetch the results of the query
            rows = cursor.fetchall()
            conn.close()
            # Convert the BLOB data to JSON dicts
            json_full[entry] = [[json.loads(row[0].decode('utf-8')), json.loads(row[1].decode('utf-8')), json.loads(row[2].decode('utf-8'))] for row in rows]

    logger.debug(2)
    logger.debug(json_full)
    return json_full


with app.app_context():
    db.create_all()


# middleware to log requests
@app.before_request
def log_request_info():
    if request.endpoint is not None and request.endpoint != 'static':
        ip_address = request.remote_addr
        username = session.get('username', 'Unknown')  # Standardwert 'Unknown', falls kein Benutzer angemeldet ist
        logger.info(f"Request from IP: {ip_address}, User: {username}, Endpoint: {request.path}")

# main page
@app.route('/')
def index():
    # go to welcome page
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
            # After successful login, redirect to the stored URL
            next_page = session.get('next', url_for('query_page'))  # Use a default if 'next' doesn't exist
            return redirect(next_page)
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
        logger.info(f"New user {username} created.")
        return redirect(url_for('login'))

    return render_template('signup.html')


# Send Query to Ontodocker
@app.route('/queryexec', methods=['POST'])
def execute_sparql_query():
    if 'username' in session:
        # Send the query to the Ontodocker and return the result
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
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('query_page')
        return redirect(url_for('login'))


# Database page (query page simple version)
@app.route('/database')
def query_page_simple():
    if 'username' in session:
        return render_template('queryPageSimple.html')
    else:
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('query_page_simple')
        return redirect(url_for('login'))


@app.route('/database_new', methods=['GET', 'POST'])
def database():
    if 'username' in session:
        conn = get_db_connection() 
        cursor = conn.cursor()

        # Check if the request method is POST
        if request.method == 'POST':
            # Get the search term from the form
            search_term = request.form.get('nameInput')

            # Modify the query to include the search term
            query = "SELECT Unique_ID, filename, type FROM uploads WHERE deleted_by_user = 0 and mapped = 1 and filename LIKE ? ORDER BY filename"
            cursor.execute(query, (f"%{search_term}%",)) # Execute the query
        else:
            query = "SELECT Unique_ID, filename, type FROM uploads WHERE deleted_by_user = 0 and mapped = 1 ORDER BY filename"
            cursor.execute(query,) # Execute the query

        # Fetch the results of the query
        rows = cursor.fetchall()

        # Convert rows to dictionaries and format uploadDate
        allData = []

        for row in rows:
            filename_without_extension = os.path.splitext(row[1])[0]
            result = {
                'Unique_ID': row[0],
                'filename': filename_without_extension,
                'type': row[2]
            }

            allData.append(result)

        # Close the connection
        conn.close()
   
        if request.method == 'POST':
            # If it's arequest, return the data in JSON format
            return jsonify(allData)
        else:
            # If it's not a request, render the template
            return render_template('database.html', data=allData)
    else:
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('database')
        return redirect(url_for('login'))



# Plotting page (query page simple version)
@app.route('/plot')
def plot_page():
    if 'username' in session:
        return render_template('plot.html')
    else:
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('plot_page')
        return redirect(url_for('login'))


# Upload page
@app.route('/upload')
def upload_page():
    if 'username' in session:
        return render_template('uploadForm.html')
    else:
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('upload_page')
        return redirect(url_for('login'))


# My Files page
@app.route('/files')
def my_files():
    if 'username' in session:
        user = session['username'] 
        conn = get_db_connection() 
        cursor = conn.cursor()

        # admin sees all files
        if user == 'admin':
            query = "SELECT Unique_ID, filename, type, uploadDate, mapped, deleted_by_user FROM uploads ORDER BY filename"
            cursor.execute(query) # Execute the query
        else:
            # user sees only his files
            query = "SELECT Unique_ID, filename, type, uploadDate, mapped, deleted_by_user FROM uploads WHERE user = ? and deleted_by_user = 0 ORDER BY filename"
            cursor.execute(query, (user,)) # Execute the query
       
        # Fetch the results of the query
        rows = cursor.fetchall()

        # Convert rows to dictionaries and format uploadDate
        uploads = []

        unmapped_files = False

        for row in rows:
            filename = row[1]

            upload = {
                'Unique_ID': row[0],
                'filename': filename,
                'type': row[2],
                'uploadDate': row[3],
                'mapped':row[4],
                'deletedByUser':row[5]
            }

            if upload['mapped'] == 0:
                unmapped_files = True

            date_object = datetime.strptime(upload['uploadDate'], '%Y-%m-%dT%H:%M:%S.%f')
            upload['uploadDate'] = date_object.strftime('%d/%m/%y %H:%M')
            uploads.append(upload)


        # Close the connection
        conn.close()

        return render_template('myFiles.html', uploads=uploads, unmapped_files=unmapped_files)
    else:
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('my_files')
        return redirect(url_for('login'))


# Admin page
@app.route('/admin')
def admin_page():
    if session['username'] == 'admin':
        return render_template('admin.html')
    else:
        return redirect(url_for('login'))


# Logout
@app.route('/logout')
def logout():
    # delete username and url out of session
    session.pop('username', None)
    session.pop('next', None)

    return redirect(url_for('index'))



@app.route('/updateDeletedByUser', methods=['POST'])
def update_deleted_by_user():
    data = request.json
    unique_id = data.get('removeFile')
    if unique_id:
        conn = get_db_connection()
        cursor = conn.cursor()
        cursor.execute("UPDATE uploads SET deleted_by_user = 1 WHERE Unique_ID = ?", (unique_id,))
        conn.commit()
        conn.close()
        return jsonify({"message": "Your delete request has been sent.","status":200})
    else:
        return jsonify({"error": "Invalid request","status":400}), 400



@app.route('/adminData', methods=['POST', 'GET'])
def get_admin_data():
    # if user is not admin, return unauthorized
    if session['username'] != 'admin':
        return jsonify({"error": "Unauthorized"}), 401

    # GET-Request
    if request.method == 'GET':
        # return config without 'DOCKER_TOKEN' key
        filtered_config = {key: value for key, value in config.items() if key != 'DOCKER_TOKEN' and key != 'SECRET_KEY'}
        users = User.query.all()
        usernames = [user.username for user in users]
        return jsonify({"config": filtered_config, "users": usernames})
    
    # POST-Request
    elif request.method == 'POST':
        data = request.json

        # clear data
        if data.get("clearData"):
            if data["clearData"]:
                clear_dataset(config)
                logger.info("Ontodocker cleared")
                try:
                    os.remove("upload.db")
                    logger.info("db cleared")
                    init_db()
                    logger.info("db initialized")
                except Exception as e:
                    logger.error(f'Error {e}')
                return jsonify({"message": "Data cleared"})
            
        # remove file (unique_id)
        elif data.get("removeFile"):
            if data["removeFile"]:
                logger.info(f"Removing file {data['removeFile']}")
                conn = get_db_connection()
                cursor = conn.cursor()

                # check if file exists in uploads
                query = "SELECT * FROM uploads WHERE Unique_ID = ?"
                cursor.execute(query, (data["removeFile"],))
                rowdata = cursor.fetchone()

                # check if file exists in uidlookup
                query = "SELECT * FROM uidlookup WHERE Unique_ID = ?"
                cursor.execute(query, (data["removeFile"],))
                rowdata2 = cursor.fetchone()

                # check if file exists in Ontodocker
                query = f'''
                SELECT ?g WHERE {{
                    ?g ?p "{data["removeFile"]}".
                }}'''
                # send query
                results = send_sparql_query(query, config)

                def check_if_key_has_empty_value(d, key):
                    """
                    Checks if a key in a potentially nested dictionary has an empty value.
                
                    :param d: The dictionary to check
                    :param key: The key to check
                    :return: True if the key has an empty value, False otherwise
                    """
                
                    for k, v in d.items():
                        if k == key:
                            return not bool(v)
                        elif isinstance(v, dict):
                            if check_if_key_has_empty_value(v, key):
                                return True
                    return False
                
                # if results is not empty, delete from Ontodocker
                if not check_if_key_has_empty_value(results, "g"):
                    if rowdata:
                        delete_specific_triples_from_endpoint(rowdata["ttl"], config)
                        logger.info(f"File {data['removeFile']} deleted from Ontodocker")
    
                # remove from both tables
                if rowdata:
                    cursor.execute("DELETE FROM uploads WHERE Unique_ID = ?", (data["removeFile"],))
                    conn.commit()
                    logger.info(f"File {data['removeFile']} deleted from uploads")
                if rowdata2:
                    cursor.execute("DELETE FROM uidlookup WHERE Unique_ID = ?", (data["removeFile"],))
                    conn.commit()
                    logger.info(f"File {data['removeFile']} deleted from uidlookup")
                conn.close()

                # if file not found
                if not rowdata and not rowdata2:
                    logger.error(f"File {data['removeFile']} not found")
                    return jsonify({"error": "File not found"}), 404
                return jsonify({"message": "File removed"})
            
        # update config
        else:
            logger.info("Updating config")
            config['ontodocker_url'] = data["config"]["ontodocker_url"]
            config["DOCKER_TOKEN"] = data["config"]["DOCKER_TOKEN"]
            config["triplestore_server"] = data["config"]["triplestore_server"]
            config["dataset_name"] = data["config"]["dataset_name"]

            with open('config.json', 'w') as file:
                json.dump(config, file, indent=4)
            logger.info("Config updated")

            for entry in data["users"]:
                # check if user exists
                user = User.query.filter_by(username=entry).first()

                if user:
                    # if user exists, delete user
                    db.session.delete(user)
                    db.session.commit()
                    logger.info(f"User {entry} deleted")
                else:
                    # if user does not exist, return error
                    return jsonify({'error': 'User not found'}), 404
                
    # if not GET or POST
    else:
        return jsonify({"error": "Method Not Allowed"}), 405


# query mixture (in upload page)
@app.route('/search-mixture', methods=['POST'])
def search_mixture():
    if 'username' in session:
        data = request.json
        uid = data['mixtureName']
        # Search, if the mixture exists in the database
        conn = get_db_connection()
        cursor = conn.cursor()

        # Sql query to search for the mixture
        query = "SELECT * FROM uidlookup WHERE Name = ? OR Unique_ID = ?"
        cursor.execute(query, (uid, uid))

        # Fetch the results of the query
        rowdata = cursor.fetchone()
        conn.close()

        if rowdata:
            return jsonify({'message': f'Mischung {uid} erfolgreich gefunden!', 'mixtureID': rowdata['Unique_ID']})
        else:
            # If the mixture does not exist in the database
            return jsonify({'message': 'Mischung nicht gefunden.'})
        
@app.route('/get-mixtures', methods=['GET'])
def get_mixtures():
    if 'username' in session:
        conn = get_db_connection()
        cursor = conn.cursor()

        # Sql query to search for the mixture
        query = "SELECT filename,Mixture_ID FROM uploads where mapped = 1 and deleted_by_user = 0 and type = 'Mixture';"
        cursor.execute(query,)

        # Fetch all the results of the query
        rowdata = cursor.fetchall()
        conn.close()

        # Extract the mixture names and ids from the row data
        mixtures = [{'name': row['filename'], 'id': row['Mixture_ID']} for row in rowdata]

        # Return the mixtures as a JSON response
        return jsonify({'mixtures': mixtures})

    else:
        return jsonify({'error': 'Nicht angemeldet'}), 403


def get_db_connection():
    try:
        conn = sqlite3.connect(upload_db)
        conn.row_factory = sqlite3.Row
    except sqlite3.Error as e:
        logger.error(f"Database connection error: {e}")
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
                url TEXT,
                blob BLOB,
                Json BLOB,
                ttl BLOB,
                Json_Specimen BLOB,
                ttl_Specimen BLOB,
                UploadDate TEXT,
                Mapped INTEGER,
                deleted_by_user BOOLEAN DEFAULT 0 NOT NULL,
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
             'Specimen': '../cpto/Specimen_KG_Template.ttl',
             'CompressiveStrength': '../cpto/CompressiveStrength_KG_Template.ttl'}

    def add_data(rowname, data):
        conn = get_db_connection()
        cursor = conn.cursor()
        # update data
        query = f"UPDATE uploads SET {rowname} = ? WHERE Unique_ID = ?"
        try:
            cursor.execute(query, (data, unique_id))
            conn.commit()
            logger.debug(f"{rowname}-Wert für {unique_id} erfolgreich aktualisiert.")
            success = 1

        # Error handling
        except Exception as ex:
            logger.error(f"Fehler beim Aktualisieren des {rowname}-Werts für {unique_id}: {ex}")
            conn.rollback()
            success = 0
        finally:
            conn.close()

        return success

    def get_data(uid=None):
        conn = get_db_connection()
        cursor = conn.cursor()

        if uid is None:
            uid = unique_id

        # sql query to get data
        query = "SELECT * FROM uploads WHERE Unique_ID = ?"
        cursor.execute(query, (uid,))

        # fetch the results of the query
        rowdata = cursor.fetchone()
        conn.close()

        return rowdata

    def upload_to_docker(data):
        success = upload_binary_to_existing_docker(data, config)
        if success != 0:
            add_data('Error', 1)

    # start of the processing
    logger.debug(f"Starte die Verarbeitung für {unique_id}")

    # fetch data
    row = get_data()

    # check if data is already mapped
    if row['Mapped'] != 0:
        logger.error(f'Error already mapped with {unique_id}')
        return

    # if json was uploaded
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
            add_data('Json', json.dumps(json_data).encode('utf-8'))
        elif row['type'] == 'EModule':
            if row['filetype'] == 'xml':
                mix_data = get_data(row['Mixture_ID'])
                json_data = xml_to_json(row['blob'], mix_data['Json'])
                emodule_json = json_data[0]
                specimen_json = json_data[1]
                emodule_json['ID'] = row['unique_id']
                specimen_json['ID'] = row['unique_id']
                emodule_json['SpecimenID'] = row['unique_id']
                specimen_json['MixtureID'] = row['Mixture_ID']
                add_data('Json', json.dumps(emodule_json).encode('utf-8'))
                add_data('Json_Specimen', json.dumps(specimen_json).encode('utf-8'))
        elif row['type'] == 'CompressiveStrength':
            if row['filetype'] == 'dat':
                # Process the raw data
                processed_data = processed_rawdata(row['blob'])
                print(processed_data)

                # Convert the processed data to a CSV 
                #csv_data = processed_data.to_csv(index=False)

    # Set ID and mixtureID in json!
    # Mapping
    # Aktualisierten Daten abfragen
    row = get_data()

    if row['Json'] == '':
        logger.error(f"No Json for {unique_id}")
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
            logger.error(f'Error in Mixture mapping: {e}')
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
            logger.error(f'Error in Placeholderreplacement: {e}')
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
                logger.error(f'Error in Placeholderreplacement: {e}')
                add_data('Error', 1)
                return

    logger.debug(f"Verarbeitung für {unique_id} abgeschlossen.")
    add_data('Mapped', 1)


# handle uploaded data
@app.route('/dataUpload', methods=['POST'])
def data_upload():
    init_db()
    file_types = ['xlsx', 'xls', 'csv', 'dat', 'txt', 'json', 'xml']
    if 'username' not in session:
        return jsonify({'error': 'Nicht angemeldet'}), 403

    if 'type' not in request.form:
        return jsonify({'error': 'Kein Typ angegeben'}), 400

    if 'Mixture_ID' not in request.form:
        return jsonify({'error': 'Keine Mixture ID angegeben'}), 400

    type = request.form['type']
    mixtureID = str(request.form['Mixture_ID'])
    user = session['username']  # get username

    # connection to the database
    conn = get_db_connection()
    cursor = conn.cursor()

    # current time
    uploaddate = datetime.now().isoformat()


    file_keys = [key for key in request.files.keys() if key.startswith('file')]
    if file_keys:  # Check if there are any files
        for file_key in file_keys:
            file = request.files[file_key]
            if file.filename != '':
                file_name = file.filename

                # Check if the file already exists in the database
                cursor.execute('SELECT * FROM uploads WHERE filename = ? and deleted_by_user = 0', (file_name,))
                data = cursor.fetchone()

                # If data is not None, then the file exist in the database
                if data is not None:
                    conn.close()
                    return jsonify({'message': "File " + file_name + " already exists! ",
                                    'status': 409}), 200
                else:
                    # extract file extension
                    _, file_extension = os.path.splitext(file.filename)
                    filetype = file_extension.lstrip('.')
                    if filetype not in file_types:
                        return jsonify({'error': 'Unsupported Type'}), 400
                    # save as blob
                    file_blob = file.read()

                    if type == 'Mixture':
                        unique_id = mixtureID
                    else:
                        unique_id = str(uuid.uuid4())

                    # insert data into the database
                    cursor.execute('INSERT INTO uploads (user, filetype, filename, type, blob, Mixture_ID, Unique_ID, UploadDate,Mapped, Error) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                                (user, filetype, file_name, type, file_blob, mixtureID, unique_id, uploaddate, 0, 0))
                    conn.commit()
                    cursor.execute('INSERT INTO uidlookup (Unique_ID, Name) VALUES (?, ?)',
                                (unique_id, file_name.split(".")[0]))
                    conn.commit()

                    # start async function for each file
                    thread = threading.Thread(target=async_function, args=(unique_id,))
                    thread.start()

    elif 'url' in request.form and request.form['url'] != '':
        url = request.form['url']
        response = requests.get(url)

        # Check if the response is HTML
        if 'text/html' in response.headers['Content-Type']:
            return jsonify({'message': "Invalid file type. Please provide a direct link to a xls, csv, txt or dat file.",
                            'status': 400}), 200
    
        file = FileStorage(BytesIO(response.content), filename=url.split('/')[-1])
        file_name = url.split('/')[-1]
        file_blob = file.read()

        # extract file extension
        _, file_extension = os.path.splitext(file.filename)
        filetype = file_extension.lstrip('.')
        if filetype not in file_types:
            return jsonify({'message': "Invalid file type. Please provide a direct link to a xls, csv, txt or dat file.",
                            'status': 400}), 200
        
        # Check if the file already exists in the database
        cursor.execute('SELECT * FROM uploads WHERE filename = ? and deleted_by_user = 0', (file_name,))
        data = cursor.fetchone()

        # If data is not None, then the file exists in the database
        if data is not None:
            conn.close()
            return jsonify({'message': "This file already exists.",
                            'status': 409}), 200
        else:
            if type == 'Mixture':
                unique_id = mixtureID
            else:
                unique_id = str(uuid.uuid4())
            # insert data into the database
            cursor.execute('INSERT INTO uploads (user, filetype, filename, url, type, blob, Mixture_ID, Unique_ID, UploadDate, '
                        'Mapped, Error) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)',
                        (user, filetype, file_name, url, type, file_blob, mixtureID, unique_id, uploaddate, 0, 0))
            conn.commit()
            cursor.execute('INSERT INTO uidlookup (Unique_ID, Name) VALUES (?, ?)',
                        (unique_id, file_name.split(".")[0]))
            conn.commit()

            # start async function
            thread = threading.Thread(target=async_function, args=(unique_id,))
            thread.start()

    else:
        return jsonify({'error': 'No Data found'}), 400

    conn.close()

    

    return jsonify({'message': "Your files have been uploaded successfully.",
                    'uniqueID': unique_id,
                    'status': 200}), 200


@app.route('/getJson', methods=['GET'])
def get_json():
    if 'username' in session:
        json_list = create_json_file()
        logger.debug(1)
        logger.debug(json_list)
        return jsonify({'json': json_list})
    else:
        return jsonify({'error': 'Nicht angemeldet'}), 403

@app.route('/rawdownload')
def raw_download():

    if 'username' not in session:
        return jsonify({'error': 'Nicht angemeldet'}), 403
    
    temp_directory = "temp/"

    if not os.path.exists(temp_directory):
        os.makedirs(temp_directory)

    if not os.path.exists("zip/"):
        os.makedirs("zip/")

    if not len(os.listdir(temp_directory)) == 0:
        shutil.rmtree(temp_directory)
        os.makedirs(temp_directory)
    
    if not len(os.listdir("zip/")) == 0:
        shutil.rmtree("zip/")
        os.makedirs("zip/")

    def zip_directory(folder_path, output_zip_path):
        # create a ZipFile object
        with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # iterate over all the files in directory
            for foldername, subfolders, filenames in os.walk(folder_path):
                for filename in filenames:
                    # create complete filepath of file in directory
                    file_path = os.path.join(foldername, filename)
                    # create complete filepath of file in zip
                    arcname = os.path.relpath(file_path, start=folder_path)
                    # Add file to zip
                    zipf.write(file_path, arcname=arcname)

        logger.info(f"ZIP-Datei wurde erstellt: {output_zip_path}")

    def get_data(uid, row="Unique_ID", table="uploads"):
        conn = get_db_connection()
        cursor = conn.cursor()

        # sql query to get data
        query = f"SELECT * FROM {table} WHERE {row} = ?"
        cursor.execute(query, (uid,))

        # fetch the results of the query
        rowdata = cursor.fetchone()
        conn.close()

        return rowdata

    # Get the file ID from the URL parameters
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
            # Path where the extracted files will be saved
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}.{row['filetype']}")

            # write BLOB data to a file
            with open(output_path, 'wb') as file:
                file.write(row['blob'])

        if row["Json"]:
            # path where the extracted files will be saved
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}.json")

            # write BLOB data to a file
            with open(output_path, 'wb') as file:
                file.write(row['Json'])

        if row["ttl"]:
            # path where the extracted files will be saved
            output_path = os.path.join(temp_directory, f"{row['filename'].split('.')[0]}.ttl")

            # write BLOB data to a file
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
            # send the zip file as a response
            return send_file(zip_output_path, as_attachment=True)
        except Exception as e:
            logger.warning(f"Error in Zip: {e}")
            # if file not found
            abort(404, description="File not found.")
    else:
        # if no id found
        abort(400, description="No file ID provided.")

if __name__ == '__main__':
    app.run(debug=True)