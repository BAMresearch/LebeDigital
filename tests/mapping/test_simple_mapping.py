from lebedigital.mapping.simple_mapping import placeholderreplacement

def test_simple_mapping ():
    # defining paths : Simple ONTOLOGY Mapped
    mappedPath = "../usecases/SimpleOntology/Mapping/SimpleOntologyMappedProbe1.ttl"  # for current ontology

    # defining paths : ONTOLOGY
    ontoPath = "../usecases/SimpleOntology/SimpleOntology.ttl"

    # defining paths : METADATA
    dataPath = "../usecases/SimpleOntology/Mapping/probe1.yaml"  # for now only working with this one, more later

    ## loading simple mapping turtle file
    with open(mappedPath, 'r') as file:
        lines = file.readlines()

    exact_result = placeholderreplacement(ontoPath, dataPath)
    assert (set(lines).discard('\n')) == (set(exact_result).discard('\n'))
