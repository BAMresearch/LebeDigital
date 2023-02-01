from lebedigital.SimpleOntology.Mapping.simple_mapping import placeholderreplacement

def test_simple_mapping ():
    print ("S T A R T - T E S T I N G")

    # defining paths : Simple ONTOLOGY Mapped
    mappedPath = "../../lebedigital/SimpleOntology/Mapping/SimpleOntologyMappedProbe1.ttl"  # for current ontology

    # defining paths : ONTOLOGY
    ontoPath = "../../lebedigital/SimpleOntology/SimpleOntology.ttl"

    # defining paths : METADATA
    dataPath = "../../lebedigital/SimpleOntology/Mapping/probe1.yaml"  # for now only working with this one, more later

    ## loading simple mapping turtle file
    with open(mappedPath, 'r') as file:
        lines = file.readlines()

    exact_result = placeholderreplacement(ontoPath, dataPath)
    assert (set(lines).discard('\n')) == (set(exact_result).discard('\n'))



test_simple_mapping()


