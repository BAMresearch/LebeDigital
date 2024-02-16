import sqlite3


# Verbindung zur Datenbank herstellen
conn = sqlite3.connect('../../jsons.db')
c = conn.cursor()

# Eine SQL-Abfrage, die alle Daten aus einer Tabelle abruft
tabelle_name = 'EModule'  # Ersetzen Sie 'meine_tabelle' mit dem tatsächlichen Tabellennamen
sql_abfrage = f'SELECT * FROM {tabelle_name}'

# Die Abfrage ausführen
c.execute(sql_abfrage)

# Alle Zeilen der Abfrageergebnisse abrufen
ergebnisse = c.fetchall()

# Ergebnisse ausgeben
for zeile in ergebnisse:
    print(zeile)

# Verbindung schließen
conn.close()
