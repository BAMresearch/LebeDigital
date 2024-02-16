import sqlite3
import json


def db_interface(daten, name):
    # Verbindung zur Datenbank herstellen
    conn = sqlite3.connect('jsons.db')
    c = conn.cursor()

    # Schlüssel direkt aus dem Dictionary extrahieren, um Spaltennamen zu bestimmen
    schluessel = daten.keys()
    spalten = ', '.join(f"{k} TEXT" for k in schluessel)

    # Tabelle dynamisch erstellen
    c.execute(f"CREATE TABLE IF NOT EXISTS {name} ({spalten})")

    # Daten einfügen
    # Da 'daten' ein einzelnes Dictionary ist, fügen wir es direkt ein
    spalten = ', '.join(daten.keys())
    platzhalter = ', '.join('?' * len(daten))
    werte = tuple(json.dumps(value) if isinstance(value, (dict, list)) else value for value in daten.values())
    c.execute(f"INSERT INTO {name} ({spalten}) VALUES ({platzhalter})", werte)

    # Änderungen commiten und Verbindung schließen
    conn.commit()
    conn.close()


def upload_to_db():
    ttl_files = {
        'Mixture': 'files/json/mix.json',
        'EModule': 'files/json/emodule.json',
        'Specimen': 'files/json/specimen.json'
    }

    for name, path in ttl_files.items():
        try:
            # Versuch, die JSON-Datei zu öffnen
            with open(path, 'r') as file:
                daten = json.load(file)
                db_interface(daten, name)
        except FileNotFoundError:
            # Datei wurde nicht gefunden
            print("Die Datei 'daten.json' wurde nicht gefunden. Überspringe das Einlesen.")
        #except Exception as e:
            # Andere Fehler
        #    print(f"Ein unerwarteter Fehler ist aufgetreten: {e}")