import rdflib
import pandas as pd
from rdflib.namespace import RDF, RDFS, OWL
from rdflib import Literal


def ont_to_excel():
    # RDF-Datei laden
    g = rdflib.Graph()
    g.parse("ConcreteOntology.owl", format="turtle")

    # Abfrage für die benötigten Informationen
    query = """
    SELECT DISTINCT ?entity ?label_en ?label_de ?def_en ?def_de ?source_en ?source_de ?remarks
    WHERE {
      ?entity a owl:Class .
      OPTIONAL {?entity rdfs:label ?label_en . FILTER (lang(?label_en) = "en")}
      OPTIONAL {?entity rdfs:label ?label_de . FILTER (lang(?label_de) = "de")}
      OPTIONAL {?entity <http://www.w3.org/ns/prov#definition> ?def_en . FILTER (lang(?def_en) = "en")}
      OPTIONAL {?entity <http://www.w3.org/ns/prov#definition> ?def_de . FILTER (lang(?def_de) = "de")}
      OPTIONAL {?entity <https://w3id.org/pmd/co/definitionSource> ?source_en . FILTER (lang(?source_en) = "en")}
      OPTIONAL {?entity <https://w3id.org/pmd/co/definitionSource> ?source_de . FILTER (lang(?source_de) = "de")}
      OPTIONAL {?entity rdfs:comment ?remarks}
    }
    """

    # Abfrageergebnisse
    results = g.query(query)

    # Daten in ein DataFrame konvertieren
    data = []
    for row in results:
        data.append([str(e) if e is not None else "nan" for e in row])

    df = pd.DataFrame(data, columns=["Entity", "Label(en)", "Label(de)", "Definition(en)", "Definition(de)", "Source(en)", "Source(de)", "Bemerkungen"])

    # Entfernen des "https://w3id.org/cpto/" Präfixes, falls vorhanden
    df['Entity'] = df['Entity'].apply(lambda x: x.replace("https://w3id.org/cpto/", "") if x.startswith("https://w3id.org/cpto/") else x)

    # Export in Excel
    df.to_excel("output.xlsx", index=False)


def excel_to_ont():
    # Excel-Tabelle einlesen
    df = pd.read_excel("output.xlsx")

    # RDF-Graph initialisieren
    g = rdflib.Graph()

    # Namespace definieren (ersetzen Sie "http://example.org/" mit Ihrem Namespace)
    ns = rdflib.Namespace("https://w3id.org/cpto")
    ns2 = rdflib.Namespace("http://www.w3.org/ns/prov#")
    pmd = rdflib.Namespace("https://w3id.org/pmd/co#")

    # Einträge aus der Excel-Tabelle in RDF umwandeln
    for _, row in df.iterrows():
        class_uri = ns[row['Entity']]

        # Eine Klasse erstellen
        g.add((class_uri, RDF.type, OWL.Class))

        # Überprüfen und Hinzufügen von Labels und Definitionen, wenn vorhanden
        if pd.notna(row['Label(en)']):
            g.add((class_uri, RDFS.label, Literal(row['Label(en)'], lang="en")))
        if pd.notna(row['Label(de)']):
            g.add((class_uri, RDFS.label, Literal(row['Label(de)'], lang="de")))
        if pd.notna(row['Definition(en)']):
            g.add((class_uri, ns2.definition, Literal(row['Definition(en)'], lang="en")))
        if pd.notna(row['Definition(de)']):
            g.add((class_uri, ns2.definition, Literal(row['Definition(de)'], lang="de")))

        # Überprüfen und Hinzufügen von Quellen und Bemerkungen
        if pd.notna(row['Source(en)']):
            g.add((class_uri, pmd.definitionSource, Literal(row['Source(en)'], lang="en")))
        if pd.notna(row['Source(de)']):
            g.add((class_uri, pmd.definitionSource, Literal(row['Source(de)'], lang="de")))
        if pd.notna(row['Bemerkungen']):
            g.add((class_uri, RDFS.comment, Literal(row['Bemerkungen'], lang="en")))

    # RDF-Datei speichern
    g.serialize(destination="updated_ontology.owl", format="turtle")


print("This tool allows to export all created Classes to the Taxonomie Table automatically.")
print('For it to work, please put the Ontology in the same folder with the name: "ConcreteOntology.owl".')
print('For Exporting select (2). The script will create the Taxonomie.xlsx as "output.xlsx"')
if input("Ont to Excel (2) or ! Experimental ! Excel to Ont (1)") == "1":
    excel_to_ont()
else:
    ont_to_excel()