import pandas as pd
import matplotlib.pyplot as plt

def visualisiere_daten(csv_dateiname, plot_type):

    # Lese den Datensatz
    df = pd.read_csv(csv_dateiname, sep=';')


    # Erhalte die Namen der ersten beiden Spalten
    spalte1, spalte2 = df.columns[:2]

    # Erstelle einen Plot
    plt.figure(figsize=(10,6))

    if plot_type == 'plot':
        plt.plot(df[spalte1], df[spalte2])
    else:
        plt.bar(df[spalte1], df[spalte2])
    # Füge Titel und Beschriftungen hinzu
    plt.title('Visualisierung von ' + spalte1 + ' gegen ' + spalte2)
    plt.xlabel(spalte1)
    plt.ylabel(spalte2)

    # Zeige den Plot
    plt.show()

# Setze hier den Namen deiner CSV-Datei ein
visualisiere_daten('Mappe1.csv', 'plot')
