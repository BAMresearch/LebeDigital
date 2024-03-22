import sqlite3
import time


def db_row_to_dict(row):
    # Verbinde zur SQLite-Datenbank
    conn = sqlite3.connect('web/upload.db')
    # Konfiguriere den Cursor, um Row-Objekte zurückzugeben
    conn.row_factory = sqlite3.Row

    cursor = conn.cursor()

    # Führe deine Abfrage aus
    cursor.execute("SELECT * FROM uploads WHERE id = ?", (row,))  # Beispiel für eine Abfrage nach einer bestimmten ID
    row = cursor.fetchone()

    if row is not None:
        # Konvertiere die Row in ein Dictionary
        row_dict = dict(row)
    else:
        print("Keine Daten gefunden.")
        row_dict = None

    # Schließe die Datenbankverbindung
    conn.close()

    return row_dict


def process_entries(entry_ids):
    # Hier deine Logik einfügen, um mit den Einträgen zu arbeiten
    print(f"Verarbeite Einträge: {entry_ids}")

    for row in entry_ids:
        dict = db_row_to_dict(row)

        # Erstellt JSON, wenn sie noch nicht existiert (in zukunft JSON Upload)
        if dict['Json'] is None:

        if dict['Json'] is not None and dict['type'] == 'Mixture':

        if dict['Json'] is not None and dict['Json_Specimen'] is not None:

def check_and_process_entries():
    connection = sqlite3.connect('web/upload.db')
    cursor = connection.cursor()

    # SQL-Abfrage, um IDs von Einträgen mit Mapped = 0 zu finden
    cursor.execute("SELECT id FROM uploads WHERE Mapped = 0")
    entries = cursor.fetchall()

    # Extrahiere die IDs aus den gefundenen Einträgen
    entry_ids = [entry[0] for entry in entries]

    if entry_ids:
        process_entries(entry_ids)

    # Schließe die Datenbankverbindung
    connection.close()


def main():
    try:
        while True:
            check_and_process_entries()
            time.sleep(30)  # Warte 30 Sekunden vor dem nächsten Durchlauf
    except KeyboardInterrupt:
        print("Skript wurde manuell gestoppt.")


if __name__ == "__main__":
    main()
