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
        print("ID | User | FileType | Type | Blob (size)")
        print("-" * 50)
        for row in rows:
            # Die Größe des BLOBs wird für die Anzeige ermittelt, um die Ausgabe lesbar zu halten
            print(f"{row[0]} | {row[1]} | {row[2]} | {row[3]} | {len(row[4])} bytes")
    else:
        print("Keine Einträge gefunden.")

    # Datenbankverbindung schließen
    conn.close()

if __name__ == '__main__':
    show_uploads()
