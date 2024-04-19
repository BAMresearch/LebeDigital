import pandas as pd

def read_excel_sheet(file_path, sheet_name):
    """
    Liest ein spezifisches Sheet aus einer Excel-Datei und gibt die Daten als DataFrame zurück.

    Parameter:
    file_path (str): Der Pfad zur Excel-Datei.
    sheet_name (str): Der Name des Sheets, das gelesen werden soll.

    Rückgabe:
    DataFrame: Die Daten des Sheets.
    """
    # Versucht, die Excel-Datei zu lesen
    try:
        data = pd.read_excel(file_path, sheet_name=sheet_name)
        print("Daten erfolgreich gelesen!")
        return data
    except Exception as e:
        print(f"Fehler beim Lesen der Datei: {e}")
        return None

# Setze hier den Pfad zur deiner Excel-Datei und den Namen des Sheets ein
file_path = '20240220_7188_M01.xls'
sheet_name = 'Rezeptur'

# Lese das Excel-Sheet
data_frame = read_excel_sheet(file_path, sheet_name)

# Ausgabe der Daten, wenn sie erfolgreich geladen wurden
if data_frame is not None:
    print(data_frame)
