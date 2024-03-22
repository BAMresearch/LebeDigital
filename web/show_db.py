import sqlite3

# Pfad zur Ihrer Datenbank
database_path = 'upload.db'

def show_uploads():
    # Verbindung zur Datenbank herstellen
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # Abfrage, um alle Einträge in der 'uploads'-Tabelle zu erhalten
    query = "SELECT * FROM uploads"
    cursor.execute(query)

    # Alle Zeilen der Abfrage abrufen
    rows = cursor.fetchall()

    # Überprüfen, ob Einträge vorhanden sind
    if rows:
        for row in rows:
            print('------------')
            for number, element in enumerate(row):
                print(number, element)
    else:
        print("Keine Einträge gefunden.")

    # Datenbankverbindung schließen
    conn.close()

if __name__ == '__main__':
    show_uploads()
