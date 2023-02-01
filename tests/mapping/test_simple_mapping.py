from pathlib import Path
import rdflib
import os
from lebedigital.SimpleOntology.Mapping.simple_mapping import placeholderreplacement

def test_simple_mapping (mappedPath, ontoPath, dataPath):
    print ("S T A R T - T E S T I N G")

    ## loading simple mapping turtle file
    with open(mappedPath, 'r') as file:
        lines = file.readlines()

    exact_result = placeholderreplacement(ontoPath, dataPath)

    assert lines == exact_result

# defining paths : ONTOLOGY
ontoPath = "../../lebedigital/SimpleOntology/SimpleOntology.ttl"

# defining paths : METADATA
dataPath = "../../lebedigital/SimpleOntology/Mapping/probe1.yaml" #for now only working with this one, more later

# defining paths : Simple ONTOLOGY Mapped
mappedPath = "../../lebedigital/SimpleOntology/Mapping/SimpleOntologyMapped.ttl" #for current ontology

test_simple_mapping(mappedPath, ontoPath, dataPath)


