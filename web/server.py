import os, sys

script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_directory, '..'))  # Add the parent directory to the path

import threading, json, uuid, sqlite3, time
from werkzeug.security import generate_password_hash, check_password_hash
from flask_sqlalchemy import SQLAlchemy
from scripts.upload.upload_script import send_sparql_query, upload_ttl_to_fuseki, clear_dataset, get_backup
from scripts.mapping.mixmapping import mappingmixture
from scripts.mapping.unit_conversion import unit_conversion_json
from scripts.mapping.mappingscript import placeholderreplacement
from scripts.rawdataextraction.emodul_xml_to_json import xml_to_json
from scripts.rawdataextraction.ComSt_xml_to_json import comSt_xml_to_json
from scripts.rawdataextraction.ComSt_generate_processed_data import processed_rawdata
from scripts.rawdataextraction.ComSt_metadata_extraction import extract_metadata_ComSt
from scripts.rawdataextraction.mixture_xml_to_json import mix_metadata
from datetime import timedelta, datetime
from flask import Flask, request, render_template, redirect, url_for, session, jsonify, abort, send_file
from loguru import logger
from pathlib import Path
import requests
from werkzeug.datastructures import FileStorage
from io import BytesIO
import zipfile
import shutil
from urllib.parse import urlparse, parse_qs
import base64
import mimetypes
import magic
import re
from collections import deque

# --------------------------------- Setup --------------------------------- #


# --- Setup Flask --- #

def load_config():
    # Check if config.json exists, if not create it from template
    if not os.path.exists('config.json'):
        with open('config.template.json', 'r') as template_file:
            config = json.load(template_file)
            config['SECRET_KEY'] = base64.b64encode(os.urandom(24)).decode('utf-8')
            with open('config.json', 'w') as config_file:
                json.dump(config, config_file, indent=4)
    else:
        with open('config.json', 'r') as config_file:
            config = json.load(config_file)
            if config['SECRET_KEY'] == 'Your Secret Key (Gets created first run)':
                config['SECRET_KEY'] = base64.b64encode(os.urandom(24)).decode('utf-8')
                with open('config.json', 'w') as file:
                    json.dump(config, file, indent=4)
    
    return config

# Initialize Flask
app = Flask(__name__)

# Load configuration
config = load_config()

# Configure Flask app
app.config['SECRET_KEY'] = config['SECRET_KEY']
app.config['APPLICATION_ROOT'] = config['APPLICATION_ROOT']
app.config['SQLALCHEMY_DATABASE_URI'] = config['DATABASE']['URI']
app.config['SESSION_PERMANENT'] = config['SESSION']['PERMANENT']
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=config['SESSION']['LIFETIME_DAYS'])
app.config['SESSION_TYPE'] = config['SESSION']['TYPE']


# --- Setup Logging --- #

# Path to the logs directory
baseDir = Path(__file__).parent
logDir = os.path.join(baseDir, "logs")

# make sure the logs directory exists
os.makedirs(logDir, exist_ok=True)

# Path to the log files
debugLogPath = os.path.join(logDir, "debug_{time}.log")
infoLogPath = os.path.join(logDir, "info_{time}.log")
errorLogPath = os.path.join(logDir, "error_{time}.log")

# Configure the logger, delete 10 days old log
logger.configure(handlers=[
    {"sink": errorLogPath, "level": "ERROR", "rotation": "5 MB", "retention": "10 days"},
    {"sink": debugLogPath, "level": "DEBUG", "rotation": "5 MB", "retention": "10 days"},
    {"sink": infoLogPath, "level": "INFO", "rotation": "5 MB", "retention": "10 days"}
])


# --- Setup User db --- #

db = SQLAlchemy(app)

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


# -- Setup Main db -- #

# Main Database
upload_db = 'main.db'

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
                additional_file BLOB,
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

init_db()


# --------------------------------- Functions --------------------------------- #

# extracts all json information from the database in two lists
def create_json_file():
    '''
    This function extracts all json information from the database.
    The first list contains the json information of the uploads table.
    The second list contains the json_specimen information of the uploads table.
    
    :return: One Json Object with keys for the types (Mixture eg) and two lists as values containing lists of all info of json and json_specimen.
    '''

    typeliste = ['Mixture', 'EModule', 'CompressiveStrength']
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
            json_full[entry] = [json.loads(row['Json'].decode('utf-8')) for row in rows if row['Json']]

        else:
            query = "SELECT u1.Json, u1.Json_specimen, u2.Json FROM uploads u1 INNER JOIN uploads u2 ON u1.Mixture_ID = u2.Unique_ID WHERE u1.type = ?"
            cursor.execute(query, (entry,))
            # Fetch the results of the query
            rows = cursor.fetchall()
            conn.close()
            # Convert the BLOB data to JSON dicts
            json_full[entry] = [[json.loads(row[0].decode('utf-8')), json.loads(row[1].decode('utf-8')), json.loads(row[2].decode('utf-8'))] for row in rows]

    logger.debug("--- Created Json Files ---")
    return json_full


# Mapping function
def async_function(unique_id):

    # Paths to the templates
    paths = {'EModule': '../cpto/EModuleOntology_KG_Template.ttl',
             'Specimen': '../cpto/Specimen_KG_Template.ttl',
             'CompressiveStrength': '../cpto/CompressiveStrength_KG_Template.ttl'}

    # ----- Functions ----- #

    # add data to the database
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

    # gets the data from the database
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

    # update the ttl in the dataset
    def upload_to_docker(data):
        success = upload_ttl_to_fuseki(data, config)
        if success != True:
            add_data('Error', 1)


    # ------ Processing ------ #

    logger.debug(f"Started processing data for: {unique_id}")

    # fetch data
    row = get_data()

    # check if data is already mapped
    if row['Mapped'] != 0:
        logger.error(f'Error already mapped with {unique_id}')
        return


    # -- Extracting the data to json -- #

    # check what filetype was uploaded (json or other)
    if row['filetype'] == 'json':
        if row['Json'] is None:
            json_data = json.loads(row['blob'].decode('utf-8'))
            json_data['ID'] = unique_id
            if row['type'] != 'Mixture':
                json_data['mixtureID'] = row['Mixture_ID']
            add_data('Json', json.dumps(json_data).encode('utf-8'))

    # if its not a json file
    else:
        # mixture xml extraction
        if row['type'] == 'Mixture':
            json_data = mix_metadata(row['blob'], row['filename'])
            json_data['ID'] = row['unique_id']
            add_data('Json', json.dumps(json_data).encode('utf-8'))
        # emodule xml extraction
        elif row['type'] == 'EModule':
            mix_data = get_data(row['Mixture_ID'])
            json_data = xml_to_json(row['blob'], mix_data['Json'])
            emodule_json = json_data[0]
            specimen_json = json_data[1]
            emodule_json['ID'] = row['unique_id']
            specimen_json['ID'] = row['unique_id']
            emodule_json['SpecimenID'] = row['unique_id']
            specimen_json['MixtureID'] = row['Mixture_ID']
            emodule = unit_conversion_json(emodule_json)
            add_data('Json', json.dumps(emodule).encode('utf-8'))
            specimen = unit_conversion_json(specimen_json)
            add_data('Json_Specimen', json.dumps(specimen).encode('utf-8'))
        # compressive strength extraction
        elif row['type'] == 'CompressiveStrength':
            mix_data = get_data(row['Mixture_ID'])
            if row['filetype'] == 'dat':
                # Process the raw data
                processed_data = processed_rawdata(row['blob'], row['unique_id'])

                # metadata extraction
                comSt_data= extract_metadata_ComSt(row['blob'], mix_data['Json'], processed_data)

            elif row['filetype'] == 'xml':
                print('entered')
                comSt_data = comSt_xml_to_json(row['blob'], mix_data['Json'])
                
            comSt_json = comSt_data[0]
            specimen_json = comSt_data[1]
            comSt_json['ID'] = row['unique_id']
            specimen_json['ID'] = row['unique_id']
            comSt_json['SpecimenID'] = row['unique_id']
            specimen_json['MixtureID'] = row['Mixture_ID']
            comSt = unit_conversion_json(comSt_json)
            add_data('Json', json.dumps(comSt).encode('utf-8'))
            specimen = unit_conversion_json(specimen_json)
            add_data('Json_Specimen', json.dumps(specimen).encode('utf-8'))


    # ---- Mapping the json to ttl ---- #

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


# --------------------------------- Routes --------------------------------- #

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
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        user = User.query.filter_by(username=username).first()
        if user and user.check_password(password):  # check password
            session['username'] = user.username
            # set session for user
            # After successful login, redirect to the stored URL
            next_page = session.get('next', url_for('database'))  # Use a default if 'next' doesn't exist
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
        return jsonify(results)


# Query page
@app.route('/query')
def query_page():
    if 'username' in session:
        return render_template('queryPage.html')
    else:
        # Store the URL the user was trying to access in the session data
        session['next'] = url_for('query_page')
        return redirect(url_for('login'))


@app.route('/database', methods=['GET', 'POST'])
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



# Plotting page
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
            query = "SELECT Unique_ID, filename, type, user, uploadDate, mapped, deleted_by_user FROM uploads ORDER BY filename"
            cursor.execute(query) # Execute the query
        else:
            # user sees only his files
            query = "SELECT Unique_ID, filename, type, user, uploadDate, mapped, deleted_by_user FROM uploads WHERE user = ? and deleted_by_user = 0 ORDER BY filename"
            cursor.execute(query, (user,)) # Execute the query
       
        # Fetch the results of the query
        rows = cursor.fetchall()

        # Convert rows to dictionaries and format uploadDate
        uploads = []

        unmapped_files = False

        for row in rows:
            upload = {
                'Unique_ID': row[0],
                'filename': row[1],
                'type': row[2],
                'uploadedBy': row[3],
                'uploadDate': row[4],
                'mapped': row[5],
                'deletedByUser': row[6]
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
        # return config without sensitive info
        filtered_config = config.copy()  # Create a copy of the config

        # Remove sensitive keys
        sensitive_keys = ['DOCKER_TOKEN', 'SECRET_KEY']
        for key in sensitive_keys:
            filtered_config.pop(key, None)

        # Remove password from fuseki config if it exists
        if 'fuseki' in filtered_config:
            fuseki_config = filtered_config['fuseki'].copy()  # Create a copy of fuseki config
            fuseki_config.pop('password', None)  # Remove password
            filtered_config['fuseki'] = fuseki_config
        
        users = User.query.all()
        usernames = [user.username for user in users]

        # Get log files
        log_data = get_log_files()
        
        return jsonify({
            "config": filtered_config, 
            "users": usernames,
            "logs": log_data
        })
    
    # POST-Request
    elif request.method == 'POST':
        data = request.json

        # Clear logs
        if data.get("clearLogs"):
            try:
                baseDir = Path(__file__).parent
                logDir = os.path.join(baseDir, "logs")
                
                # Get all log files
                log_files = [f for f in os.listdir(logDir) if f.endswith('.log')]
                
                # Delete each log file
                for file in log_files:
                    file_path = os.path.join(logDir, file)
                    try:
                        os.remove(file_path)
                    except Exception as e:
                        logger.error(f"Error deleting log file {file}: {str(e)}")
                
                # Create new log files (this will happen automatically when logging resumes)
                logger.info("All log files cleared by admin")
                
                return jsonify({"message": "All logs cleared successfully"})
            except Exception as e:
                error_msg = f"Error clearing logs: {str(e)}"
                logger.error(error_msg)
                return jsonify({"error": error_msg}), 500


        # clear data
        if data.get("clearData"):
            if data["clearData"]:
                clear_dataset(config)
                logger.info("Fuseki cleared")
                try:
                    os.remove(upload_db)
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
                
                cursor.execute("DELETE FROM uploads WHERE Unique_ID = ?", (data["removeFile"],))
                conn.commit()
                logger.info(f"File {data['removeFile']} deleted from uploads")
            
                cursor.execute("DELETE FROM uidlookup WHERE Unique_ID = ?", (data["removeFile"],))
                conn.commit()
                logger.info(f"File {data['removeFile']} deleted from uidlookup")
                conn.close()

                return jsonify({"message": "The file has been deleted permanently."})
            
        # update config
        else:
            logger.info("Updating config")
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


def get_log_files():
    """Get the most recent log files and their content"""
    baseDir = Path(__file__).parent
    logDir = os.path.join(baseDir, "logs")
    log_data = {}
    
    try:
        # List all files in the logs directory
        files = os.listdir(logDir)
        
        # Filter for .log files and sort by modification time (most recent first)
        log_files = [f for f in files if f.endswith('.log')]
        log_files.sort(key=lambda x: os.path.getmtime(os.path.join(logDir, x)), reverse=True)
        
        # Get the content of each log file (limit to most recent files)
        for file in log_files:  # Limit to 5 most recent files
            file_path = os.path.join(logDir, file)
            file_size = os.path.getsize(file_path)
            
            # Skip empty files
            if file_size == 0:
                continue

            try:
                with open(file_path, 'r') as f:
                    # Read last 1000 lines (or adjust as needed)
                    content = list(deque(f, 1000))
                    log_data[file] = {
                        'content': content,
                        'size': os.path.getsize(file_path),
                        'modified': os.path.getmtime(file_path)
                    }
            except Exception as e:
                log_data[file] = {
                    'error': f"Could not read file: {str(e)}",
                    'size': os.path.getsize(file_path),
                    'modified': os.path.getmtime(file_path)
                }
                
    except Exception as e:
        return {'error': f"Could not access log directory: {str(e)}"}
    
    return log_data


@app.route('/backup', methods=['GET'])
def backup():
    # Authentifizierung prüfen
    if 'username' not in session or session['username'] != 'admin':
        return jsonify({'error': 'Nicht angemeldet'}), 403

    # Backup der externen Quelle abrufen
    backup_stream = get_backup(config)
    if not backup_stream:
        return jsonify({'error': 'Backup von der externen Quelle fehlgeschlagen.'}), 500

    

    try:
        # Dateien in eine ZIP-Datei bündeln
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
            # Lokale Datei hinzufügen
            zip_file.write("main.db", "main.db")
            # Externes Backup hinzufügen
            zip_file.writestr("external_backup.nt", backup_stream.getvalue())

        # Pufferinhalt als ZIP-Datei zurückgeben
        zip_buffer.seek(0)
        return send_file(
            zip_buffer,
            mimetype='application/zip',
            as_attachment=True,
            download_name='backup_files.zip'
        )
    except Exception as e:
        logger.error(f"Fehler beim Erstellen der ZIP-Datei: {str(e)}")
        return jsonify({"error": str(e)}), 500



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


@app.route('/database_individual', methods=['POST'])
def get_individual_data():
    if 'username' in session:
        unique_id = request.json.get('id')

        conn = get_db_connection()
        cursor = conn.cursor()

        # Sql query to get the data of the unique_id
        query = "SELECT Json, Json_Specimen FROM uploads WHERE Unique_ID = ?"
        cursor.execute(query, (unique_id,))

        # Fetch the results of the query
        rowdata = cursor.fetchone()
        # Initialize empty dictionaries for JSON data
        json_data = {}
        json_specimen_data = {}

        # Check if rowdata is not None
        if rowdata:
            # Load Json data if it exists
            if rowdata[0]:
                json_data = json.loads(rowdata[0])
            
            # Load Json_Specimen data if it exists
            if rowdata[1]:
                json_specimen_data = json.loads(rowdata[1])

        # Merge the two JSON objects
        merged_json = {**json_data, **json_specimen_data}

        conn.close()

        # Return the data as a JSON response
        return jsonify(merged_json)

    else:
        return jsonify({'error': 'Nicht angemeldet'}), 403
    

def download_file_from_url(url):
    parsed_url = urlparse(url)
    domain = parsed_url.netloc

    # Handle specific cases for GitHub, Dropbox, SharePoint, Google Drive, etc.
    if 'view.officeapps.live.com' in domain:
        query_params = parse_qs(parsed_url.query)
        if 'src' in query_params:
            url = query_params['src'][0]
        else:
            return None, None, None, "Invalid URL. Please provide a valid link."
    elif domain == 'github.com':
        url = url.replace('github.com', 'raw.githubusercontent.com').replace('/blob', '')
    elif 'dropbox.com' in domain:
        if 'dl=0' in url:
            url = url.replace('dl=0', 'dl=1')
        elif '?dl=0' not in url and '?raw=1' not in url:
            url += '?dl=1'
        else:
            return None, None, None, "Invalid Dropbox URL. Please provide a valid link."
    elif 'drive.google.com' in domain or 'docs.google.com' in domain:   
        if 'edit' in url:
            url = url.replace('edit', 'export')
        elif 'edit' not in url:
            url += '/export'
        else:
            return None, None, None, "Invalid Google Drive URL. Please provide a valid link."
    else:
        # Handle other URLs
        return None, None, None, "Sorry! Please provide a valid link from GitHub, Dropbox or Google Drive."
    
    response = requests.get(url, allow_redirects=True)
    
    content_type = response.headers.get('Content-Type', '')
    if 'text/html' in content_type:
        return None, None, None, "Sorry! We could not access your file. Please download it and upload the file from your device."

    # Correctly unpack the filename and extension returned from `get_filename_from_response`
    file_name, extension = get_filename_from_response(response, url)
    
    # Create the file storage object
    file = FileStorage(BytesIO(response.content), filename=file_name)
    
    return file, file_name, extension, None


def clean_filename(file_name):
    # First, split the filename and extension
    base_name, extension = os.path.splitext(file_name)

    # If the base_name itself contains multiple dots, keep only the first part
    if '.' in base_name:
        base_name = base_name.split('.')[0]

    # Remove any leading dot from the extension
    extension = extension.lstrip('.')
    
    # Reassemble the cleaned filename with its valid extension
    cleaned_filename = f"{base_name}.{extension}" if extension else base_name
    
    return cleaned_filename, extension


def get_filename_from_response(response, url):
    # Try to get the filename from the Content-Disposition header
    content_disposition = response.headers.get('Content-Disposition')
    if content_disposition:
        # Check for both standard and extended filename formats
        fname = re.findall(r'filename\*?=(?:UTF-8\'\')?(.+?)(?:;|$)', content_disposition)
        if len(fname) > 0:
            # Decode filename if it's encoded (e.g., UTF-8)
            file_name = fname[0].strip().strip('"').strip("'")

            # Clean the filename to remove extra parts after the extension
            return clean_filename(file_name)
    
    # If not found, use the last part of the URL
    file_name = url.split('/')[-1].split('?')[0]
    
    # Clean the filename to remove extra parts after the extension
    return clean_filename(file_name)

# handle uploaded data
@app.route('/dataUpload', methods=['POST'])
def data_upload():
    init_db()
    file_types = ['dat', 'json', 'xml', 'xls', 'xlsx']
    if 'username' not in session:
        return jsonify({'error': 'Please Login first.'}), 403

    if 'type' not in request.form:
        return jsonify({'error': 'Please select a type.'}), 400

    if 'Mixture_ID' not in request.form:
        return jsonify({'error': 'No Mixture ID found!'}), 400

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
                if check_file_exists(file_name, cursor):
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

                    # Call async_function directly instead of in a thread
                    async_function(unique_id)

                    # start async function for each file
                    #thread = threading.Thread(target=async_function, args=(unique_id,))
                    #thread.start()

    elif 'url' in request.form and request.form['url'] != '':
        url = request.form['url']
        file, file_name, extension, error_message = download_file_from_url(url)

        if error_message:
            return jsonify({'message': error_message, 'status': 400}), 200

        file_blob = file.read()
        file.seek(0)  # Reset file pointer to the beginning

        filetype = extension
        if filetype not in file_types:
            return jsonify({'message': f"Invalid file type: {filetype}. Please provide a xls, csv, txt or dat file.",
                            'status': 400}), 200
        
        # Check if the file already exists in the database
        if check_file_exists(file_name, cursor):
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

            # Call async_function directly instead of in a thread
            async_function(unique_id)

            # start async function for each file
            #thread = threading.Thread(target=async_function, args=(unique_id,))
            #thread.start()
    else:
        return jsonify({'error': 'No Data found'}), 400
    
    conn.close()
    
    return jsonify({'message': "Your files have been uploaded successfully.",
                    'uniqueID': unique_id,
                    'status': 200}), 200


#manually upload mixture
@app.route('/new_mixture')
def new_mixture():
    return render_template('newMixture.html')


@app.route('/submit_mixture', methods=['POST'])
def submit_mixture():
    user = session.get('username')  # get username
    if not user:
        return jsonify({"status": 403, "message": "User not authenticated"}), 403
    
    try:
        # Access all form data
        form_data = request.form.to_dict()
        
        # Filter out keys with empty values
        filtered_form_data = {key: value for key, value in form_data.items() if value}

        # Additional fixed values
        additional = {
            "Binder1_Density_Unit": "kg/dm^3",
            "Binder1_Amount_Unit": "kg/m^3",
            "Water1_Density_Unit": "kg/dm^3",
            "Water1_Total_Unit": "kg/m^3",
            "Water1_Effective_Unit": "kg/m^3",
            "Aggregate1_Density_Unit": "kg/dm^3",
            "Aggregate1_Amount_Unit": "kg/m^3",
            "Addition1_Density_Unit": "kg/dm^3",
            "Addition1_Amount_Unit": "kg/m^3",
            "Admixture1_Density_Unit": "kg/dm^3",
            "Admixture1_Amount_Unit": "kg/m^3",
            "Fiber1_Density_Unit": "kg/dm^3",
            "Fiber1_Amount_Unit": "kg/m^3",
            "Air1_Density_Unit": "kg/dm^3",
            "Air1_Amount_Unit": "kg/m^3",
            "RawDataFile": "Download"
        }

        # Update the filtered_form_data with unit values
        filtered_form_data.update(additional)

        # current time
        upload_date = datetime.now().isoformat()
        type = "Mixture" 
        mixtureID = str(uuid.uuid4())
        filename = f"{filtered_form_data['HumanReadableID']}.json"
        fileType = 'json'
        print(filename)

        # sort the JSON keys
        def custom_sort_key(key):
            # Prioritize Lab, MixingDate, and HumanReadableID
            priority_keys = ["Lab", "MixingDate", "HumanReadableID", "RawDataFile"]
            if key in priority_keys:
                return (priority_keys.index(key), key)
            # Custom sorting for Binder, Water, etc.
            elif key.startswith("Binder"):
                return (len(priority_keys) + 1, key)
            elif key.startswith("Water"):
                return (len(priority_keys) + 2, key)
            elif key.startswith("Aggregate"):
                return (len(priority_keys) + 3, key)
            elif key.startswith("Addition"):
                return (len(priority_keys) + 4, key)
            elif key.startswith("Admixture"):
                return (len(priority_keys) + 5, key)
            elif key.startswith("Fiber"):
                return (len(priority_keys) + 6, key)
            elif key.startswith("Air"):
                return (len(priority_keys) + 7, key)
            else:
                return (len(priority_keys) + 8, key)  # Other keys come last

        mix_dict = dict(sorted(filtered_form_data.items(), key=lambda item: custom_sort_key(item[0])))
        mix_dict = unit_conversion_json(mix_dict)

        # correcting wrong spelled keys

        if mix_dict.get('Water1_Ratio'):
            mix_dict['WaterCementRatio'] = mix_dict['Water1_Ratio']
            mix_dict.pop('Water1_Ratio')

        if mix_dict.get('HumanReadableID'):
            mix_dict['humanreadableID'] = mix_dict['HumanReadableID']
            mix_dict.pop('HumanReadableID')

        mix_dict['ID'] = mixtureID

        json_blob = json.dumps(mix_dict).encode('utf-8')
        
        # Handle file upload and save as BLOB in the database
        additional_file_blob = None
        file = request.files.get('file')
        
        if file:
            additional_file_blob = file.read()  # Read file content as binary
            file_info = {
                'filename': file.filename,
                'content': base64.b64encode(additional_file_blob).decode('utf-8')

            }
            additional_file_blob = json.dumps(file_info).encode('utf-8')


        # Connect to the database
        try:
            conn = get_db_connection()
            cursor = conn.cursor()

            # Check if the filename already exists
            if check_file_exists(filename, cursor):
                conn.close()
                return jsonify({'message': f"Human-readable ID {filtered_form_data['HumanReadableID']} already exists!", 'status': 409}), 409
        
            cursor.execute('''
                INSERT INTO uploads 
                (user, filetype, filename, type, blob, Json, Mixture_ID, Unique_ID, additional_file, UploadDate, Mapped, Error) 
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (user, fileType, filename, type, json_blob, json_blob, mixtureID, mixtureID, additional_file_blob,  upload_date, 0, 0))
            conn.commit()
            # start async function
            async_function(mixtureID)
            # thread = threading.Thread(target=async_function, args=(mixtureID,))
            # thread.start()
        except Exception as e:
            conn.rollback()
            return jsonify({"status": 500, "message": f"Database error: {str(e)}"}), 500
        finally:
            conn.close()
        
        # Example response
        return jsonify({"success": True, "data": form_data})
        
    except Exception as e:
        return jsonify({"success": False, "message": str(e)})


#manually upload compressiveStrength
@app.route('/new_compressive_strength')
def new_compressive_strength():
    return render_template('newCompressiveStrength.html')

@app.route('/submit_compressive_strength', methods=['POST'])
def submit_compressive_strength():
    user = session.get('username')  # get username
    if not user:
        return jsonify({"status": 403, "message": "User not authenticated"}), 403
    
    # Retrieve the file and JSON data from the form
    file = request.files.get('file')
    comst_json = request.form.get('comst')
    specimen_json = request.form.get('specimen')

    # Check if JSON data is provided
    if not comst_json or not specimen_json:
        return jsonify({"status": 400, "message": "No JSON data provided"}), 400
    
    comst = json.loads(comst_json)
    specimen = json.loads(specimen_json)

    # current time
    upload_date = datetime.now().isoformat()

    UniqueID = str(uuid.uuid4())
    comst['ID'] = comst['specimenID'] = UniqueID
    specimen['ID'] = UniqueID

    filename = f"{comst['humanreadableID']}.json"
    type = 'CompressiveStrength'

    comst = unit_conversion_json(comst)
    specimen = unit_conversion_json(specimen)

    # Convert jsons to a JSON string and then to bytes
    comst_json = json.dumps(comst).encode('utf-8')
    specimen_json = json.dumps(specimen).encode('utf-8')

    # Handle file upload and save as BLOB in the database
    additional_file_blob = None
    if file:
        additional_file_blob = file.read()  # Read file content as binary
        file_info = {
            'filename': file.filename,
            'content': base64.b64encode(additional_file_blob).decode('utf-8')

        }
        additional_file_blob = json.dumps(file_info).encode('utf-8')


    # Connect to the database
    try:
        conn = get_db_connection()
        cursor = conn.cursor()

        # Check if the filename already exists
        if check_file_exists(filename, cursor):
            conn.close()
            return jsonify({'message': f"Human-readable ID {comst['humanreadableID']} already exists!", 'status': 409}), 409
        
        cursor.execute('''
            INSERT INTO uploads 
            (user, filetype, filename, type, blob, Json, Json_Specimen, additional_file, Mixture_ID, Unique_ID, UploadDate, Mapped, Error) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (user, "json", filename, type, comst_json, comst_json, specimen_json, additional_file_blob, specimen['MixtureID'], UniqueID, upload_date, 0, 0))
        conn.commit()
        async_function(UniqueID)
        # start async function
        #thread = threading.Thread(target=async_function, args=(UniqueID,))
        #thread.start()
    except Exception as e:
        logger.debug(e)
        conn.rollback()
        return jsonify({"status": 500, "message": f"Database error: {str(e)}"}), 500
    finally:
        conn.close()
    
    return jsonify({"success": True})


#manually upload emodule
@app.route('/new_emodule')
def new_emodule():
    return render_template('newEmodule.html')

@app.route('/submit_emodule', methods=['POST'])
def submit_emodule():
    user = session.get('username')  # get username
    if not user:
        return jsonify({"status": 403, "message": "User not authenticated"}), 403
    
    # Retrieve the file and JSON data from the form
    file = request.files.get('file')
    emodule_json = request.form.get('emodule')
    specimen_json = request.form.get('specimen')

    # Check if JSON data is provided
    if not emodule_json or not specimen_json:
        return jsonify({"status": 400, "message": "No JSON data provided"}), 400
    
    emodule = json.loads(emodule_json)
    specimen = json.loads(specimen_json)

    # current time
    upload_date = datetime.now().isoformat()

    UniqueID = str(uuid.uuid4())
    emodule['ID'] = emodule['specimenID'] = UniqueID
    specimen['ID'] = UniqueID

    filename = f"{emodule['humanreadableID']}.json"
    type = 'EModule'

    emodule = unit_conversion_json(emodule)
    specimen = unit_conversion_json(specimen)

    # Convert jsons to a JSON string and then to bytes
    emodule_json = json.dumps(emodule).encode('utf-8')
    specimen_json = json.dumps(specimen).encode('utf-8')

    # Handle file upload and save as BLOB in the database
    additional_file_blob = None
    if file:
        additional_file_blob = file.read()  # Read file content as binary
        file_info = {
            'filename': file.filename,
            'content': base64.b64encode(additional_file_blob).decode('utf-8')

        }
        additional_file_blob = json.dumps(file_info).encode('utf-8')


    # Connect to the database
    try:
        conn = get_db_connection()
        cursor = conn.cursor()

        # Check if the filename already exists
        if check_file_exists(filename, cursor):
            conn.close()
            return jsonify({'message': f"Human-readable ID {emodule['humanreadableID']} already exists!", 'status': 409}), 409
        
        cursor.execute('''
            INSERT INTO uploads 
            (user, filetype, filename, type, blob, Json, Json_Specimen, additional_file, Mixture_ID, Unique_ID, UploadDate, Mapped, Error) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (user, "json", filename, type, emodule_json, emodule_json, specimen_json, additional_file_blob, specimen['MixtureID'], UniqueID, upload_date, 0, 0))
        conn.commit()
        
        async_function(UniqueID)

        # start async function
        #thread = threading.Thread(target=async_function, args=(UniqueID,))
        #thread.start()
    except Exception as e:
        logger.debug(e)
        conn.rollback()
        return jsonify({"status": 500, "message": f"Database error: {str(e)}"}), 500
    finally:
        conn.close()
    
    return jsonify({"success": True})


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

        if row["additional_file"]:
            file_info = json.loads(row["additional_file"].decode('utf-8'))
            file_content = base64.b64decode(file_info['content'])
            original_filename = file_info['filename']

            # Detect file type
            file_type = magic.from_buffer(file_content, mime=True)
            
            # Ensure the filename has the correct extension
            file_extension = mimetypes.guess_extension(file_type)
            if file_extension and not original_filename.lower().endswith(file_extension.lower()):
                original_filename += file_extension

            # Path where the additional file will be saved
            output_path = os.path.join(temp_directory, original_filename)

            # Write additional_file data to a file
            with open(output_path, 'wb') as file:
                file.write(file_content)

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


def check_file_exists(filename, cursor):
    """
    Check if the file already exists in the database.

    Args:
    - filename (str): The name of the file to check.
    - cursor (sqlite3.Cursor): The database cursor to execute the query.

    Returns:
    - bool: True if the file exists, False otherwise.
    """
    cursor.execute('SELECT * FROM uploads WHERE filename = ? AND deleted_by_user = 0', (filename,))
    data = cursor.fetchone()
    return data is not None


# ---------------------------- #

if __name__ == '__main__':
    app.run(debug=True)
