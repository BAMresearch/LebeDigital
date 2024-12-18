import re
from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.namespace import RDF
import json
from scripts.mapping.unit_conversion import unit_conversion
from loguru import logger

'''
This Script is for the mapping of only the mixture, since it is fundamentally different from everything else.
In the mixture new classes and values need to be dynamically mapped and created. With the placeholder method
this is not really possible. That's why this script creates the ttl file from the ground up, meaning that all changes
need to be done here. This is not perfect but it works just fine. This Script is made so changes shouldn't relly
be necessary and changes to the main Ontology shouldn't create the need to change this Script, only renaming of classes.
If Classes will be renamed just Search and Replace them ( ! Upper and Lowercase ! )
If the structure of the KG gets changed slightly you can change it here, its sorted.
If the structure fundamentally changes this needs to be changed here too.
'''


def mappingmixture(blob_data):
    """
    Mappes the give data to the Mixture KG

    :param blob_data: json file in binary format containing all relevant informations
    :return: ttl file as binary
    """

    # Versuche, den BLOB als JSON zu laden
    try:
        # Der BLOB ist binär, konvertiere ihn zuerst in einen String
        json_data = blob_data.decode('utf-8')
        # Lade den JSON-Inhalt
        metadata = json.loads(json_data)
        # Verwende deine Funktion `unit_conversion` auf die geladenen Metadaten
        metadata = unit_conversion(metadata)
    except Exception as e:
        print(f'Error reading Mixture json: {e}')
        return

    # Create a new graph
    g = Graph()

    # Define the namespaces
    co = Namespace("https://w3id.org/pmd/co/")
    cpto = Namespace("https://w3id.org/cpto/")
    qudt = Namespace("http://qudt.org/schema/qudt/")
    ns4 = Namespace("http://purl.org/spar/datacite/")
    owl = Namespace("http://www.w3.org/2002/07/owl#")
    rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
    xsd = Namespace("http://www.w3.org/2001/XMLSchema#")

    # bind prefix
    g.bind("co", co)
    g.bind("cpto", cpto)
    g.bind("qudt", qudt)
    g.bind("ns4", ns4)
    g.bind("owl", owl)
    g.bind("rdfs", rdfs)
    g.bind("xsd", xsd)

    # add the basic triples, no logic here
    g.add((ns4.hasIdentifier, RDF.type, owl.ObjectProperty))
    g.add((xsd.anyURI, RDF.type, rdfs.Datatype))
    g.add((xsd.dateTime, RDF.type, rdfs.Datatype))
    g.add((xsd.decimal, RDF.type, rdfs.Datatype))
    g.add((xsd.float, RDF.type, rdfs.Datatype))
    g.add((xsd.string, RDF.type, rdfs.Datatype))
    g.add((cpto[""], RDF.type, owl.Ontology))
    g.add((co.characteristic, RDF.type, owl.ObjectProperty))
    g.add((co.composedOf, RDF.type, owl.ObjectProperty))
    g.add((co.input, RDF.type, owl.ObjectProperty))
    g.add((co.unit, RDF.type, owl.ObjectProperty))
    g.add((co.value, RDF.type, owl.ObjectProperty))
    g.add((cpto.Addition, RDF.type, owl.Class))
    g.add((cpto.Admixture, RDF.type, owl.Class))
    g.add((cpto.AggregateSize, RDF.type, owl.Class))
    g.add((cpto.Cement, RDF.type, owl.Class))
    g.add((cpto.Air, RDF.type, owl.Class))
    g.add((cpto.Water, RDF.type, owl.Class))
    g.add((cpto.WaterCementRatio, RDF.type, owl.Class))
    g.add((cpto.MaterialComposition, RDF.type, owl.Class))
    g.add((cpto.Content, RDF.type, owl.Class))
    g.add((cpto.Manufacturer, RDF.type, owl.Class))
    g.add((cpto.Kind, RDF.type, owl.Class))
    g.add((cpto.Date, RDF.type, owl.Class))
    g.add((cpto.Product, RDF.type, owl.Class))
    g.add((cpto.RelativeDensity, RDF.type, owl.Class))
    g.add((cpto.Aggregate, RDF.type, owl.Class))
    g.add((co.Dataset, RDF.type, owl.Class))
    g.add((co.Laboratory, RDF.type, owl.Class))
    g.add((co.NodeInfo, RDF.type, owl.Class))
    g.add((co.Time, RDF.type, owl.Class))
    g.add((co.ProvidedIdentifier, RDF.type, owl.Class))
    g.add((co.BaseMaterial, RDF.type, owl.Class))
    g.add((qudt.Unit, RDF.type, owl.Class))

    # add the basic logic of the concrete mixture
    concreteMixture = URIRef(cpto + f"ConcreteMixture_{metadata['ID']}")

    # define the concrete mixture individual
    g.add((concreteMixture, RDF.type, owl.NamedIndividual))
    g.add((concreteMixture, RDF.type, cpto.MaterialComposition))

    # add properties
    g.add((concreteMixture, ns4.hasIdentifier, URIRef(cpto + f"ID_{metadata['ID']}")))
    g.add((concreteMixture, ns4.hasIdentifier, URIRef(cpto + f"humanreadableID_{metadata['ID']}")))
    g.add((concreteMixture, co.characteristic, URIRef(cpto + f"ExperimentInfo_{metadata['ID']}")))
    g.add((concreteMixture, co.characteristic, URIRef(cpto + f"MaterialComposition_{metadata['ID']}")))
    g.add((concreteMixture, co.characteristic, URIRef(cpto + f"WaterCementRatio_{metadata['ID']}")))
    g.add((concreteMixture, co.input, URIRef(cpto + f"RawDataFile_{metadata['ID']}")))

    # add all the units from the qudt, that are in the json

    # first go through all the units in the metadata
    for key, value in metadata.items():
        # if the key is a unit then add the unit to the ontology
        if "_Unit" in key:
            g.add((URIRef(f"{value}"), RDF.type, qudt.Unit))
            g.add((URIRef(f"{value}"), RDF.type, owl.NamedIndividual))

    # now add the cement, if multiple exist, add multiple

    # define the search pattern
    pattern = r"Binder(\d+)"

    # iterate through the keys, to find all matching keys, that are different (how many cements were used)
    found_numbers = {match for key in metadata.keys() for match in re.findall(pattern, key)}
    logger.debug(f'Found {len(found_numbers)} different Binders.')

    # create a list of the amount of Cements
    cements = list(found_numbers)

    # now create the triples for every cement
    for entry in cements:
        # add triple
        cement_individual = URIRef(cpto + f"Binder{entry}_{metadata['ID']}")
        cement_content = URIRef(cpto + f"Binder{entry}_Amount_{metadata['ID']}")
        cement_date = URIRef(cpto + f"Binder{entry}_Date_{metadata['ID']}")
        cement_density = URIRef(cpto + f"Binder{entry}_Density_{metadata['ID']}")
        cement_kind = URIRef(cpto + f"Binder{entry}_Kind_{metadata['ID']}")
        cement_manufacturer = URIRef(cpto + f"Binder{entry}_Manufacturer_{metadata['ID']}")
        cement_type = URIRef(cpto + f"Binder{entry}_Type_{metadata['ID']}")

        # cement
        g.add((cement_individual, RDF.type, co.BaseMaterial))
        g.add((cement_individual, RDF.type, owl.NamedIndividual))
        g.add((cement_individual, co.characteristic, cement_content))
        g.add((cement_individual, co.characteristic, cement_date))
        g.add((cement_individual, co.characteristic, cement_density))
        g.add((cement_individual, co.characteristic, cement_kind))
        g.add((cement_individual, co.characteristic, cement_manufacturer))
        g.add((cement_individual, co.composedOf, cement_type))

        # cement content
        g.add((cement_content, RDF.type, cpto.Content))
        g.add((cement_content, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Binder{entry}_Amount_Unit'):
            g.add((cement_content, co.unit, URIRef(f"{metadata[f'Binder{entry}_Amount_Unit']}")))
        if metadata.get(f'Binder{entry}_Amount'):
            g.add((cement_content, co.value, Literal(f"{metadata[f'Binder{entry}_Amount']}", datatype=xsd.float)))

        # cement date
        g.add((cement_date, RDF.type, cpto.Date))
        g.add((cement_date, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Binder{entry}_Date'):
            g.add((cement_date, co.value, Literal(f"{metadata[f'Binder{entry}_Date']}", datatype=xsd.string)))

        # cement density
        g.add((cement_density, RDF.type, cpto.RelativeDensity))
        g.add((cement_density, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Binder{entry}_Density_Unit'):
            g.add((cement_density, co.unit, URIRef(f"{metadata[f'Binder{entry}_Density_Unit']}")))
        if metadata.get(f'Binder{entry}_Density'):
            g.add((cement_density, co.value, Literal(f"{metadata[f'Binder{entry}_Density']}", datatype=xsd.float)))

        # cement kind
        g.add((cement_kind, RDF.type, cpto.Kind))
        g.add((cement_kind, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Binder{entry}_Kind'):
            g.add((cement_kind, co.value, Literal(f"{metadata[f'Binder{entry}_Kind']}", datatype=xsd.string)))

        # cement manufacturer
        g.add((cement_manufacturer, RDF.type, cpto.Manufacturer))
        g.add((cement_manufacturer, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Binder{entry}_Manufacturer'):
            g.add((cement_manufacturer, co.value, Literal(f"{metadata[f'Binder{entry}_Manufacturer']}", datatype=xsd.string)))

        # cement type
        g.add((cement_type, RDF.type, cpto.Cement))
        g.add((cement_type, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Binder{entry}_Type'):
            g.add((cement_type, co.value, Literal(f"{metadata[f'Binder{entry}_Type']}", datatype=xsd.string)))

    # now we add the aggregate, basically like the cement

    # define the search pattern
    pattern = r"Aggregate(\d+)"

    # iterate through the keys, to find all matching keys, that are different (how many cements were used)
    found_numbers = {match for key in metadata.keys() for match in re.findall(pattern, key)}
    logger.debug(f'Found {len(found_numbers)} different Aggregates.')

    # create a list of the amount of Cements
    aggregates = list(found_numbers)

    # now create the triples for every cement
    for entry in aggregates:
        # add triple
        aggregate_individual = URIRef(cpto + f"Aggregate{entry}_{metadata['ID']}")
        aggregate_content = URIRef(cpto + f"Aggregate{entry}_Amount_{metadata['ID']}")
        aggregate_date = URIRef(cpto + f"Aggregate{entry}_Date_{metadata['ID']}")
        aggregate_size = URIRef(cpto + f"Aggregate{entry}_Fraction_{metadata['ID']}")
        aggregate_kind = URIRef(cpto + f"Aggregate{entry}_Kind_{metadata['ID']}")
        aggregate_manufacturer = URIRef(cpto + f"Aggregate{entry}_Manufacturer_{metadata['ID']}")
        aggregate_type = URIRef(cpto + f"Aggregate{entry}_Type_{metadata['ID']}")
        aggregate_density = URIRef(cpto + f"Aggregate{entry}_Density_{metadata['ID']}")

        # aggregate
        g.add((aggregate_individual, RDF.type, cpto.Aggregate))
        g.add((aggregate_individual, RDF.type, owl.NamedIndividual))
        g.add((aggregate_individual, co.characteristic, aggregate_content))
        g.add((aggregate_individual, co.characteristic, aggregate_date))
        g.add((aggregate_individual, co.characteristic, aggregate_kind))
        g.add((aggregate_individual, co.characteristic, aggregate_manufacturer))
        g.add((aggregate_individual, co.characteristic, aggregate_density))
        g.add((aggregate_individual, co.characteristic, aggregate_size))
        g.add((aggregate_individual, co.composedOf, aggregate_type))

        # aggregate content
        g.add((aggregate_content, RDF.type, cpto.Content))
        g.add((aggregate_content, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Amount_Unit'):
            g.add((aggregate_content, co.unit, URIRef(f"{metadata[f'Aggregate{entry}_Amount_Unit']}")))
        if metadata.get(f'Aggregate{entry}_Amount'):
            g.add((aggregate_content, co.value, Literal(f"{metadata[f'Aggregate{entry}_Amount']}", datatype=xsd.float)))
        
        # aggregate date
        g.add((aggregate_date, RDF.type, cpto.Date))
        g.add((aggregate_date, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Date'):
            g.add((aggregate_date, co.value, Literal(f"{metadata[f'Aggregate{entry}_Date']}", datatype=xsd.string)))
        
        # aggregate size
        g.add((aggregate_size, RDF.type, cpto.AggregateSize))
        g.add((aggregate_size, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Fraction_Unit'):
            g.add((aggregate_size, co.unit, URIRef(f"{metadata[f'Aggregate{entry}_Fraction_Unit']}")))
        if metadata.get(f'Aggregate{entry}_Fraction'):
            g.add((aggregate_size, co.value, Literal(f"{metadata[f'Aggregate{entry}_Fraction']}", datatype=xsd.string)))


        # aggregate kind
        g.add((aggregate_kind, RDF.type, cpto.Kind))
        g.add((aggregate_kind, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Kind'):
            g.add((aggregate_kind, co.value, Literal(f"{metadata[f'Aggregate{entry}_Kind']}", datatype=xsd.string)))

        # aggregate manufacturer
        g.add((aggregate_manufacturer, RDF.type, cpto.Manufacturer))
        g.add((aggregate_manufacturer, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Manufacturer'):
            g.add((aggregate_manufacturer, co.value, Literal(f"{metadata[f'Aggregate{entry}_Manufacturer']}", datatype=xsd.string)))

        # aggregate type
        g.add((aggregate_type, RDF.type, cpto.Aggregate))
        g.add((aggregate_type, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Type'):
            g.add((aggregate_type, co.value, Literal(f"{metadata[f'Aggregate{entry}_Type']}", datatype=xsd.string)))


        # aggregate density
        g.add((aggregate_density, RDF.type, cpto.RelativeDensity))
        g.add((aggregate_density, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Aggregate{entry}_Density_Unit'):
            g.add((aggregate_density, co.unit, URIRef(f"{metadata[f'Aggregate{entry}_Density_Unit']}")))
        if metadata.get(f'Aggregate{entry}_Density'):
            g.add((aggregate_density, co.value, Literal(f"{metadata[f'Aggregate{entry}_Density']}", datatype=xsd.float)))


    # now admixture

    # define the search pattern
    pattern = r"Admixture(\d+)"

    # iterate through the keys, to find all matching keys, that are different (how many cements were used)
    found_numbers = {match for key in metadata.keys() for match in re.findall(pattern, key)}
    logger.debug(f'Found {len(found_numbers)} different Admixtures.')

    # create a list of the amount of Cements
    admixtures = list(found_numbers)

    # now create the triples for every cement
    for entry in admixtures:
        # add triple
        admixture_individual = URIRef(cpto + f"Admixture{entry}_{metadata['ID']}")
        admixture_content = URIRef(cpto + f"Admixture{entry}_Amount_{metadata['ID']}")
        admixture_date = URIRef(cpto + f"Admixture{entry}_Date_{metadata['ID']}")
        admixture_density = URIRef(cpto + f"Admixture{entry}_Density_{metadata['ID']}")
        admixture_kind = URIRef(cpto + f"Admixture{entry}_Kind_{metadata['ID']}")
        admixture_manufacturer = URIRef(cpto + f"Admixture{entry}_Manufacturer_{metadata['ID']}")
        admixture_productname = URIRef(cpto + f"Admixture{entry}_ProductName_{metadata['ID']}")
        admixture_type = URIRef(cpto + f"Admixture{entry}_Type_{metadata['ID']}")

        # cement
        g.add((admixture_individual, RDF.type, co.BaseMaterial))
        g.add((admixture_individual, RDF.type, owl.NamedIndividual))
        g.add((admixture_individual, co.characteristic, admixture_content))
        g.add((admixture_individual, co.characteristic, admixture_date))
        g.add((admixture_individual, co.characteristic, admixture_density))
        g.add((admixture_individual, co.characteristic, admixture_kind))
        g.add((admixture_individual, co.characteristic, admixture_manufacturer))
        g.add((admixture_individual, co.characteristic, admixture_productname))
        g.add((admixture_individual, co.composedOf, admixture_type))

        ## AB HIER
        # cement content
        g.add((admixture_content, RDF.type, cpto.Content))
        g.add((admixture_content, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_Amount_Unit'):
            g.add((admixture_content, co.unit, URIRef(f"{metadata[f'Admixture{entry}_Amount_Unit']}")))
        if metadata.get(f'Admixture{entry}_Amount'):
            g.add((admixture_content, co.value, Literal(f"{metadata[f'Admixture{entry}_Amount']}", datatype=xsd.float)))

        # cement date
        g.add((admixture_date, RDF.type, cpto.Date))
        g.add((admixture_date, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_Date'):
            g.add((admixture_date, co.value, URIRef(f"{metadata[f'Admixture{entry}_Date']}")))

        # cement density
        g.add((admixture_density, RDF.type, cpto.RelativeDensity))
        g.add((admixture_density, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_Density_Unit'):
            g.add((admixture_density, co.unit, URIRef(f"{metadata[f'Admixture{entry}_Density_Unit']}")))
        if metadata.get(f'Admixture{entry}_Density'):
            g.add((admixture_density, co.value, Literal(f"{metadata[f'Admixture{entry}_Density']}", datatype=xsd.float)))

        # cement kind
        g.add((admixture_kind, RDF.type, cpto.Kind))
        g.add((admixture_kind, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_Kind'):
            g.add((admixture_kind, co.value, URIRef(f"{metadata[f'Admixture{entry}_Kind']}")))

        # cement manufacturer
        g.add((admixture_manufacturer, RDF.type, cpto.Manufacturer))
        g.add((admixture_manufacturer, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_Manufacturer'):
            g.add((admixture_manufacturer, co.value, URIRef(f"{metadata[f'Admixture{entry}_Manufacturer']}")))

        # cement productname
        g.add((admixture_productname, RDF.type, co.ProvidedIdentifier))
        g.add((admixture_productname, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_ProductName'):
            g.add((admixture_productname, co.value, URIRef(f"{metadata[f'Admixture{entry}_ProductName']}")))

        # cement type
        g.add((admixture_type, RDF.type, cpto.Admixture))
        g.add((admixture_type, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Admixture{entry}_Type'):
            g.add((admixture_type, co.value, Literal(f"{metadata[f'Admixture{entry}_Type']}", datatype=xsd.string)))


    # now fibers

    # define the search pattern
    pattern = r"Fiber(\d+)"

    # iterate through the keys, to find all matching keys, that are different (how many cements were used)
    found_numbers = {match for key in metadata.keys() for match in re.findall(pattern, key)}
    logger.debug(f'Found {len(found_numbers)} different Fibers.')

    # create a list of the amount of fibers
    admixtures = list(found_numbers)

    # now create the triples for every cement
    for entry in admixtures:
        # add triple
        admixture_individual = URIRef(cpto + f"Fiber{entry}_{metadata['ID']}")
        admixture_content = URIRef(cpto + f"Fiber{entry}_Amount_{metadata['ID']}")
        admixture_date = URIRef(cpto + f"Fiber{entry}_Date_{metadata['ID']}")
        admixture_density = URIRef(cpto + f"Fiber{entry}_Density_{metadata['ID']}")
        admixture_kind = URIRef(cpto + f"Fiber{entry}_Kind_{metadata['ID']}")
        admixture_manufacturer = URIRef(cpto + f"Fiber{entry}_Manufacturer_{metadata['ID']}")
        admixture_productname = URIRef(cpto + f"Fiber{entry}_ProductName_{metadata['ID']}")
        admixture_type = URIRef(cpto + f"Fiber{entry}_Type_{metadata['ID']}")

        # cement
        g.add((admixture_individual, RDF.type, co.BaseMaterial))
        g.add((admixture_individual, RDF.type, owl.NamedIndividual))
        g.add((admixture_individual, co.characteristic, admixture_content))
        g.add((admixture_individual, co.characteristic, admixture_date))
        g.add((admixture_individual, co.characteristic, admixture_density))
        g.add((admixture_individual, co.characteristic, admixture_kind))
        g.add((admixture_individual, co.characteristic, admixture_manufacturer))
        g.add((admixture_individual, co.characteristic, admixture_productname))
        g.add((admixture_individual, co.composedOf, admixture_type))

        ## AB HIER
        # cement content
        g.add((admixture_content, RDF.type, cpto.Content))
        g.add((admixture_content, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Amount_Unit'):
            g.add((admixture_content, co.unit, URIRef(f"{metadata[f'Fiber{entry}_Amount_Unit']}")))
        if metadata.get(f'Fiber{entry}_Amount'):
            g.add((admixture_content, co.value, Literal(f"{metadata[f'Fiber{entry}_Amount']}", datatype=xsd.float)))

        # cement date
        g.add((admixture_date, RDF.type, cpto.Date))
        g.add((admixture_date, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Date'):
            g.add((admixture_date, co.value, URIRef(f"{metadata[f'Fiber{entry}_Date']}")))

        # cement density
        g.add((admixture_density, RDF.type, cpto.RelativeDensity))
        g.add((admixture_density, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Density_Unit'):
            g.add((admixture_density, co.unit, URIRef(f"{metadata[f'Fiber{entry}_Density_Unit']}")))
        if metadata.get(f'Fiber{entry}_Density'):
            g.add((admixture_density, co.value, Literal(f"{metadata[f'Fiber{entry}_Density']}", datatype=xsd.float)))

        # cement kind
        g.add((admixture_kind, RDF.type, cpto.Kind))
        g.add((admixture_kind, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Kind'):
            g.add((admixture_kind, co.value, URIRef(f"{metadata[f'Fiber{entry}_Kind']}")))

        # cement manufacturer
        g.add((admixture_manufacturer, RDF.type, cpto.Manufacturer))
        g.add((admixture_manufacturer, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Manufacturer'):
            g.add((admixture_manufacturer, co.value, URIRef(f"{metadata[f'Fiber{entry}_Manufacturer']}")))

        # cement productname
        g.add((admixture_productname, RDF.type, co.ProvidedIdentifier))
        g.add((admixture_productname, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Name'):
            g.add((admixture_productname, co.value, URIRef(f"{metadata[f'Fiber{entry}_Name']}")))

        # cement type
        g.add((admixture_type, RDF.type, cpto.Admixture))
        g.add((admixture_type, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Fiber{entry}_Type'):
            g.add((admixture_type, co.value, Literal(f"{metadata[f'Fiber{entry}_Type']}", datatype=xsd.string)))


    # lastly addition

    # define the search pattern
    pattern = r"Addition(\d+)"

    # iterate through the keys, to find all matching keys, that are different (how many cements were used)
    found_numbers = {match for key in metadata.keys() for match in re.findall(pattern, key)}
    logger.debug(f'Found {len(found_numbers)} different Admixtures.')

    # create a list of the amount of Cements
    additions = list(found_numbers)

    # now create the triples for every cement
    for entry in additions:
        # add triple
        addition_individual = URIRef(cpto + f"Addition{entry}_{metadata['ID']}")
        addition_content = URIRef(cpto + f"Addition{entry}_Amount_{metadata['ID']}")
        addition_date = URIRef(cpto + f"Addition{entry}_Date_{metadata['ID']}")
        addition_density = URIRef(cpto + f"Addition{entry}_Density_{metadata['ID']}")
        addition_kind = URIRef(cpto + f"Addition{entry}_Kind_{metadata['ID']}")
        addition_manufacturer = URIRef(cpto + f"Addition{entry}_Manufacturer_{metadata['ID']}")
        addition_product = URIRef(cpto + f"Addition{entry}_Product_{metadata['ID']}")
        addition_solidcontent = URIRef(cpto + f"Addition{entry}_SolidContent_{metadata['ID']}")
        addition_type = URIRef(cpto + f"Addition{entry}_Type_{metadata['ID']}")

        # cement
        g.add((addition_individual, RDF.type, co.BaseMaterial))
        g.add((addition_individual, RDF.type, owl.NamedIndividual))
        g.add((addition_individual, co.characteristic, addition_content))
        g.add((addition_individual, co.characteristic, addition_date))
        g.add((addition_individual, co.characteristic, addition_density))
        g.add((addition_individual, co.characteristic, addition_kind))
        g.add((addition_individual, co.characteristic, addition_manufacturer))
        g.add((addition_individual, co.characteristic, addition_product))
        g.add((addition_individual, co.characteristic, addition_solidcontent))
        g.add((addition_individual, co.composedOf, addition_type))

        # cement content
        g.add((addition_content, RDF.type, cpto.Content))
        g.add((addition_content, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Amount_Unit'):
            g.add((addition_content, co.unit, URIRef(f"{metadata[f'Addition{entry}_Amount_Unit']}")))
        if metadata.get(f'Addition{entry}_Amount'):
            g.add((addition_content, co.value, Literal(f"{metadata[f'Addition{entry}_Amount']}", datatype=xsd.float)))

        # cement date
        g.add((addition_date, RDF.type, cpto.Date))
        g.add((addition_date, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Date'):
            g.add((addition_date, co.value, URIRef(f"{metadata[f'Addition{entry}_Date']}")))

        # cement density
        g.add((addition_density, RDF.type, cpto.RelativeDensity))
        g.add((addition_density, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Density_Unit'):
            g.add((addition_density, co.unit, URIRef(f"{metadata[f'Addition{entry}_Density_Unit']}")))
        if metadata.get(f'Addition{entry}_Density'):
            g.add((addition_density, co.value, Literal(f"{metadata[f'Addition{entry}_Density']}", datatype=xsd.float)))

        # cement kind
        g.add((addition_kind, RDF.type, cpto.Kind))
        g.add((addition_kind, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Kind'):
            g.add((addition_kind, co.value, URIRef(f"{metadata[f'Addition{entry}_Kind']}")))

        # cement manufacturer
        g.add((addition_manufacturer, RDF.type, cpto.Manufacturer))
        g.add((addition_manufacturer, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Manufacturer'):
            g.add((addition_manufacturer, co.value, URIRef(f"{metadata[f'Addition{entry}_Manufacturer']}")))

        # cement product
        g.add((addition_product, RDF.type, cpto.Product))
        g.add((addition_product, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Product'):
            g.add((addition_product, co.value, URIRef(f"{metadata[f'Addition{entry}_Product']}")))

        # cement solidcontent
        g.add((addition_solidcontent, RDF.type, cpto.Content))
        g.add((addition_solidcontent, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_SolidContent'):
            g.add((addition_solidcontent, co.value, URIRef(f"{metadata[f'Addition{entry}_SolidContent']}")))

        # cement type
        g.add((addition_type, RDF.type, cpto.Addition))
        g.add((addition_type, RDF.type, owl.NamedIndividual))
        if metadata.get(f'Addition{entry}_Type'):
            g.add((addition_type, co.value, Literal(f"{metadata[f'Addition{entry}_Type']}", datatype=xsd.string)))

    # now water, without multiple, just one

    # add triple
    water_individual = URIRef(cpto + f"Water_{metadata['ID']}")
    water_content = URIRef(cpto + f"Water_Total_{metadata['ID']}")
    water_density = URIRef(cpto + f"Water_Density_{metadata['ID']}")
    water_type = URIRef(cpto + f"Water_Type_{metadata['ID']}")
    water_effective = URIRef(cpto + f"Water_Effective_{metadata['ID']}")

    # water
    g.add((water_individual, RDF.type, co.BaseMaterial))
    g.add((water_individual, RDF.type, owl.NamedIndividual))
    g.add((water_individual, co.characteristic, water_content))
    g.add((water_individual, co.characteristic, water_density))
    g.add((water_individual, co.characteristic, water_effective))
    g.add((water_individual, co.composedOf, water_type))

    # water content
    g.add((water_content, RDF.type, cpto.Content))
    g.add((water_content, RDF.type, owl.NamedIndividual))
    if metadata.get('Water1_Total_Unit'):
        g.add((water_content, co.unit, URIRef(f"{metadata[f'Water1_Total_Unit']}")))
    if metadata.get('Water1_Total'):
        g.add((water_content, co.value, Literal(f"{metadata[f'Water1_Total']}", datatype=xsd.float)))

    # water density
    g.add((water_density, RDF.type, cpto.RelativeDensity))
    g.add((water_density, RDF.type, owl.NamedIndividual))
    if metadata.get('Water1_Density_Unit'):
        g.add((water_density, co.unit, URIRef(f"{metadata[f'Water1_Density_Unit']}")))
    if metadata.get('Water1_Density'):
        g.add((water_density, co.value, Literal(f"{metadata[f'Water1_Density']}", datatype=xsd.float)))

    # water type
    g.add((water_type, RDF.type, cpto.Addition))
    g.add((water_type, RDF.type, owl.NamedIndividual))
    if metadata.get('Water1_Type'):
        g.add((water_type, co.value, Literal(f"{metadata[f'Water1_Type']}", datatype=xsd.string)))

    # water effective
    g.add((water_effective, RDF.type, cpto.Content))
    g.add((water_effective, RDF.type, owl.NamedIndividual))
    if metadata.get('Water1_Effective'):
        g.add((water_effective, co.value, Literal(f"{metadata[f'Water1_Effective']}", datatype=xsd.string)))

    # humanreadableID

    humanreadable_individual = URIRef(cpto + f"humanreadableID_{metadata['ID']}")

    g.add((humanreadable_individual, RDF.type, co.ProvidedIdentifier))
    g.add((humanreadable_individual, RDF.type, owl.NamedIndividual))
    g.add((humanreadable_individual, co.value, Literal(f"{metadata['humanreadableID']}", datatype=xsd.string)))

    # water cement ratio

    wz_individual = URIRef(cpto + f"WaterCementRatio_{metadata['ID']}")

    g.add((wz_individual, RDF.type, cpto.WaterCementRatio))
    g.add((wz_individual, RDF.type, owl.NamedIndividual))
    if metadata.get('WaterCementRatio'):
        g.add((wz_individual, co.value, Literal(f"{metadata['WaterCementRatio']}", datatype=xsd.float)))

    # rawdata file

    raw_individual = URIRef(cpto + f"RawDataFile_{metadata['ID']}")

    g.add((raw_individual, RDF.type, co.Dataset))
    g.add((raw_individual, RDF.type, owl.NamedIndividual))
    if metadata.get('RawDataFile'):
        g.add((raw_individual, co.value, Literal(f"{metadata['RawDataFile']}", datatype=xsd.string)))

    # mixing date

    mixing_individual = URIRef(cpto + f"MixingDate_{metadata['ID']}")

    g.add((mixing_individual, RDF.type, co.Time))
    g.add((mixing_individual, RDF.type, owl.NamedIndividual))
    if metadata.get('MixingDate'):
        g.add((mixing_individual, co.value, Literal(f"{metadata['MixingDate']}", datatype=xsd.string)))

    # lab

    lab_individual = URIRef(cpto + f"Laboratory_{metadata['ID']}")

    g.add((lab_individual, RDF.type, co.Laboratory))
    g.add((lab_individual, RDF.type, owl.NamedIndividual))
    if metadata.get('Lab'):
        g.add((lab_individual, co.value, Literal(f"{metadata['Lab']}", datatype=xsd.string)))


    # air

    air_individual = URIRef(cpto + f"Air_{metadata['ID']}")
    air_amount = URIRef(cpto + f"Air_Amount_{metadata['ID']}")
    air_density = URIRef(cpto + f"Air_Density_{metadata['ID']}")

    g.add((air_individual, RDF.type, co.Air))
    g.add((air_individual, RDF.type, owl.NamedIndividual))
    g.add((air_individual, co.characteristic, air_amount))
    g.add((air_individual, co.characteristic, air_density))

    if metadata.get('Air1_Amount'):
        g.add((air_amount, co.value, Literal(f"{metadata['Air1_Amount']}", datatype=xsd.string)))
    if metadata.get('Air1_Amount_Unit'):
        g.add((air_amount, co.unit, Literal(f"{metadata['Air1_Amount_Unit']}", datatype=xsd.string)))

    if metadata.get('Air1_Density'):
        g.add((air_density, co.value, Literal(f"{metadata['Air1_Density']}", datatype=xsd.string)))
    if metadata.get('Air1_Density_Unit'):
        g.add((air_density, co.unit, Literal(f"{metadata['Air1_Density_Unit']}", datatype=xsd.string)))

    # id

    id_individual = URIRef(cpto + f"ID_{metadata['ID']}")

    g.add((id_individual, RDF.type, co.ProvidedIdentifier))
    g.add((id_individual, RDF.type, owl.NamedIndividual))
    g.add((id_individual, co.value, Literal(f"{metadata['ID']}", datatype=xsd.string)))

    # experiment info

    experiment_individual = URIRef(cpto + f"ExperimentInfo_{metadata['ID']}")

    g.add((experiment_individual, RDF.type, co.NodeInfo))
    g.add((experiment_individual, RDF.type, owl.NamedIndividual))
    g.add((experiment_individual, co.input, lab_individual))
    g.add((experiment_individual, co.input, mixing_individual))

    # material composition (important)

    composition_individual = URIRef(cpto + f"MaterialComposition_{metadata['ID']}")

    g.add((composition_individual, RDF.type, cpto.MaterialComposition))
    g.add((composition_individual, RDF.type, owl.NamedIndividual))
    g.add((composition_individual, co.composedOf, water_individual))

    # check if multiples exist
    for entry in additions:
        g.add((composition_individual, co.composedOf, URIRef(cpto + f"Addition{entry}_{metadata['ID']}")))

    # check if multiples exist
    for entry in admixtures:
        g.add((composition_individual, co.composedOf, URIRef(cpto + f"Admixture{entry}_{metadata['ID']}")))

    # check if multiples exist
    for entry in aggregates:
        g.add((composition_individual, co.composedOf, URIRef(cpto + f"Aggregate{entry}_{metadata['ID']}")))

    # check if multiples exist
    for entry in cements:
        g.add((composition_individual, co.composedOf, URIRef(cpto + f"Binder{entry}_{metadata['ID']}")))

    try:
        # Serialisiert den Graphen im Turtle-Format und kodiert ihn in UTF-8 als binäre Daten
        serialized_data = g.serialize(format="turtle").encode('UTF-8')
        return serialized_data
    except Exception as e:
        # Fange mögliche Fehler beim Serialisieren und Kodieren ab
        print(f"Fehler beim Serialisieren des Graphen: {e}")
        return None
