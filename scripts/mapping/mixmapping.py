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


def mappingmixture(inputpath, outputpath):

    with open(inputpath, 'r') as file:
        try:
            metadata = json.load(file)
            metadata = unit_conversion(metadata)
        except Exception as e:
            print(f'Error reading Mixture json: {e}')

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
    g.add((cpto.Water, RDF.type, owl.Class))
    g.add((cpto.WaterCementRatio, RDF.type, owl.Class))
    g.add((cpto.MaterialComposition, RDF.type, owl.Class))
    g.add((cpto.Content, RDF.type, owl.Class))
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
    pattern = r"Cement(\d+)"

    # iterate through the keys, to find all matching keys, that are different (how many cements were used)
    found_numbers = {match for key in metadata.keys() for match in re.findall(pattern, key)}
    logger.debug(f'Found {len(found_numbers)} different Cements.')

    # create a list of the amount of Cements
    cements = list(found_numbers)

    # now create the triples for every cement
    for entry in cements:
        # add triple
        cement_individual = URIRef(cpto + f"Cement{entry}_{metadata['ID']}")
        cement_content = URIRef(cpto + f"Cement{entry}_Content_{metadata['ID']}")
        cement_density = URIRef(cpto + f"Cement{entry}_Density_{metadata['ID']}")
        cement_type = URIRef(cpto + f"Cement{entry}_Type_{metadata['ID']}")

        # cement
        g.add((cement_individual, RDF.type, co.BaseMaterial))
        g.add((cement_individual, RDF.type, owl.NamedIndividual))
        g.add((cement_individual, co.characteristic, cement_content))
        g.add((cement_individual, co.characteristic, cement_density))
        g.add((cement_individual, co.composedOf, cement_type))

        # cement content
        g.add((cement_content, RDF.type, cpto.Content))
        g.add((cement_content, RDF.type, owl.NamedIndividual))
        g.add((cement_content, co.unit, URIRef(f"{metadata[f'Cement{entry}_Content_Unit']}")))
        g.add((cement_content, co.value, Literal(f"{metadata[f'Cement{entry}_Content']}", datatype=xsd.float)))

        # cement density
        g.add((cement_density, RDF.type, cpto.RelativeDensity))
        g.add((cement_density, RDF.type, owl.NamedIndividual))
        g.add((cement_density, co.unit, URIRef(f"{metadata[f'Cement{entry}_Density_Unit']}")))
        g.add((cement_density, co.value, Literal(f"{metadata[f'Cement{entry}_Density']}", datatype=xsd.float)))

        # cement type
        g.add((cement_type, RDF.type, cpto.Cement))
        g.add((cement_type, RDF.type, owl.NamedIndividual))
        g.add((cement_type, co.value, Literal(f"{metadata[f'Cement{entry}_Type']}", datatype=xsd.string)))

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
        aggregate_content = URIRef(cpto + f"Aggregate{entry}_Content_{metadata['ID']}")
        aggregate_density = URIRef(cpto + f"Aggregate{entry}_Density_{metadata['ID']}")
        aggregate_size = URIRef(cpto + f"Aggregate{entry}_Size_{metadata['ID']}")
        aggregate_type = URIRef(cpto + f"Aggregate{entry}_Type_{metadata['ID']}")

        # aggregate
        g.add((aggregate_individual, RDF.type, cpto.Aggregate))
        g.add((aggregate_individual, RDF.type, owl.NamedIndividual))
        g.add((aggregate_individual, co.characteristic, aggregate_content))
        g.add((aggregate_individual, co.characteristic, aggregate_density))
        g.add((aggregate_individual, co.characteristic, aggregate_size))
        g.add((aggregate_individual, co.composedOf, aggregate_type))

        # aggregate content
        g.add((aggregate_content, RDF.type, cpto.Content))
        g.add((aggregate_content, RDF.type, owl.NamedIndividual))
        g.add((aggregate_content, co.unit, URIRef(f"{metadata[f'Aggregate{entry}_Content_Unit']}")))
        g.add((aggregate_content, co.value, Literal(f"{metadata[f'Aggregate{entry}_Content']}", datatype=xsd.float)))

        # aggregate density
        g.add((aggregate_density, RDF.type, cpto.RelativeDensity))
        g.add((aggregate_density, RDF.type, owl.NamedIndividual))
        g.add((aggregate_density, co.unit, URIRef(f"{metadata[f'Aggregate{entry}_Density_Unit']}")))
        g.add((aggregate_density, co.value, Literal(f"{metadata[f'Aggregate{entry}_Density']}", datatype=xsd.float)))

        # aggregate size
        g.add((aggregate_size, RDF.type, cpto.AggregateSize))
        g.add((aggregate_size, RDF.type, owl.NamedIndividual))
        g.add((aggregate_size, co.unit, URIRef(f"{metadata[f'Aggregate{entry}_Size_Unit']}")))
        g.add((aggregate_size, co.value, Literal(f"{metadata[f'Aggregate{entry}_Size']}", datatype=xsd.float)))

        # aggregate type
        g.add((aggregate_type, RDF.type, cpto.Aggregate))
        g.add((aggregate_type, RDF.type, owl.NamedIndividual))
        g.add((aggregate_type, co.value, Literal(f"{metadata[f'Aggregate{entry}_Type']}", datatype=xsd.string)))

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
        admixture_content = URIRef(cpto + f"Admixture{entry}_Content_{metadata['ID']}")
        admixture_density = URIRef(cpto + f"Admixture{entry}_Density_{metadata['ID']}")
        admixture_type = URIRef(cpto + f"Admixture{entry}_Type_{metadata['ID']}")

        # cement
        g.add((admixture_individual, RDF.type, co.BaseMaterial))
        g.add((admixture_individual, RDF.type, owl.NamedIndividual))
        g.add((admixture_individual, co.characteristic, admixture_content))
        g.add((admixture_individual, co.characteristic, admixture_density))
        g.add((admixture_individual, co.composedOf, admixture_type))

        # cement content
        g.add((admixture_content, RDF.type, cpto.Content))
        g.add((admixture_content, RDF.type, owl.NamedIndividual))
        g.add((admixture_content, co.unit, URIRef(f"{metadata[f'Admixture{entry}_Content_Unit']}")))
        g.add((admixture_content, co.value, Literal(f"{metadata[f'Admixture{entry}_Content']}", datatype=xsd.float)))

        # cement density
        g.add((admixture_density, RDF.type, cpto.RelativeDensity))
        g.add((admixture_density, RDF.type, owl.NamedIndividual))
        g.add((admixture_density, co.unit, URIRef(f"{metadata[f'Admixture{entry}_Density_Unit']}")))
        g.add((admixture_density, co.value, Literal(f"{metadata[f'Admixture{entry}_Density']}", datatype=xsd.float)))

        # cement type
        g.add((admixture_type, RDF.type, cpto.Admixture))
        g.add((admixture_type, RDF.type, owl.NamedIndividual))
        g.add((admixture_type, co.value, Literal(f"{metadata[f'Admixture{entry}_Type']}", datatype=xsd.string)))

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
        addition_content = URIRef(cpto + f"Addition{entry}_Content_{metadata['ID']}")
        addition_density = URIRef(cpto + f"Addition{entry}_Density_{metadata['ID']}")
        addition_type = URIRef(cpto + f"Addition{entry}_Type_{metadata['ID']}")

        # cement
        g.add((addition_individual, RDF.type, co.BaseMaterial))
        g.add((addition_individual, RDF.type, owl.NamedIndividual))
        g.add((addition_individual, co.characteristic, addition_content))
        g.add((addition_individual, co.characteristic, addition_density))
        g.add((addition_individual, co.composedOf, addition_type))

        # cement content
        g.add((addition_content, RDF.type, cpto.Content))
        g.add((addition_content, RDF.type, owl.NamedIndividual))
        g.add((addition_content, co.unit, URIRef(f"{metadata[f'Addition{entry}_Content_Unit']}")))
        g.add((addition_content, co.value, Literal(f"{metadata[f'Addition{entry}_Content']}", datatype=xsd.float)))

        # cement density
        g.add((addition_density, RDF.type, cpto.RelativeDensity))
        g.add((addition_density, RDF.type, owl.NamedIndividual))
        g.add((addition_density, co.unit, URIRef(f"{metadata[f'Addition{entry}_Density_Unit']}")))
        g.add((addition_density, co.value, Literal(f"{metadata[f'Addition{entry}_Density']}", datatype=xsd.float)))

        # cement type
        g.add((addition_type, RDF.type, cpto.Addition))
        g.add((addition_type, RDF.type, owl.NamedIndividual))
        g.add((addition_type, co.value, Literal(f"{metadata[f'Addition{entry}_Type']}", datatype=xsd.string)))

    # now water, without multiple, just one

    # add triple
    water_individual = URIRef(cpto + f"Water_{metadata['ID']}")
    water_content = URIRef(cpto + f"Water_Content_{metadata['ID']}")
    water_density = URIRef(cpto + f"Water_Density_{metadata['ID']}")
    water_type = URIRef(cpto + f"Water_Type_{metadata['ID']}")

    # cement
    g.add((water_individual, RDF.type, co.BaseMaterial))
    g.add((water_individual, RDF.type, owl.NamedIndividual))
    g.add((water_individual, co.characteristic, water_content))
    g.add((water_individual, co.characteristic, water_density))
    g.add((water_individual, co.composedOf, water_type))

    # cement content
    g.add((water_content, RDF.type, cpto.Content))
    g.add((water_content, RDF.type, owl.NamedIndividual))
    g.add((water_content, co.unit, URIRef(f"{metadata[f'Water_Content_Unit']}")))
    g.add((water_content, co.value, Literal(f"{metadata[f'Water_Content']}", datatype=xsd.float)))

    # cement density
    g.add((water_density, RDF.type, cpto.RelativeDensity))
    g.add((water_density, RDF.type, owl.NamedIndividual))
    g.add((water_density, co.unit, URIRef(f"{metadata[f'Water_Density_Unit']}")))
    g.add((water_density, co.value, Literal(f"{metadata[f'Water_Density']}", datatype=xsd.float)))

    # cement type
    g.add((water_type, RDF.type, cpto.Addition))
    g.add((water_type, RDF.type, owl.NamedIndividual))
    g.add((water_type, co.value, Literal(f"{metadata[f'Water_Type']}", datatype=xsd.string)))


    # humanreadableID

    humanreadable_individual = URIRef(cpto + f"humanreadableID_{metadata['ID']}")

    g.add((humanreadable_individual, RDF.type, co.ProvidedIdentifier))
    g.add((humanreadable_individual, RDF.type, owl.NamedIndividual))
    g.add((humanreadable_individual, co.value, Literal(f"{metadata['humanreadableID']}", datatype=xsd.string)))

    # water cement ratio

    wz_individual = URIRef(cpto + f"WaterCementRatio_{metadata['ID']}")

    g.add((wz_individual, RDF.type, cpto.WaterCementRatio))
    g.add((wz_individual, RDF.type, owl.NamedIndividual))
    g.add((wz_individual, co.value, Literal(f"{metadata['WaterCementRatio']}", datatype=xsd.decimal)))

    # rawdata file

    raw_individual = URIRef(cpto + f"RawDataFile_{metadata['ID']}")

    g.add((raw_individual, RDF.type, co.Dataset))
    g.add((raw_individual, RDF.type, owl.NamedIndividual))
    g.add((raw_individual, co.value, Literal(f"{metadata['RawDataFile']}", datatype=xsd.string)))

    # mixing date

    mixing_individual = URIRef(cpto + f"MixingDate_{metadata['ID']}")

    g.add((mixing_individual, RDF.type, co.Time))
    g.add((mixing_individual, RDF.type, owl.NamedIndividual))
    g.add((mixing_individual, co.value, Literal(f"{metadata['MixingDate']}", datatype=xsd.string)))

    # lab

    lab_individual = URIRef(cpto + f"Laboratory_{metadata['ID']}")

    g.add((lab_individual, RDF.type, co.Laboratory))
    g.add((lab_individual, RDF.type, owl.NamedIndividual))
    g.add((lab_individual, co.value, Literal(f"{metadata['Lab']}", datatype=xsd.string)))

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
        g.add((composition_individual, co.composedOf, URIRef(cpto + f"Cement{entry}_{metadata['ID']}")))

    # Speichern des Graphen in einer TTL-Datei
    with open(outputpath, "wb") as ttl_file:
        ttl_file.write(g.serialize(format="turtle").encode('UTF-8'))
        logger.debug('Mixture Mapping success.')
