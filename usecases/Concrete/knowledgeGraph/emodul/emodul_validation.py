from operator import imod
import subprocess
from time import sleep
from turtle import pos
from typing import Tuple

from torch import true_divide
import requests
import os
import sys
from pathlib import Path
from git import Repo
from rdflib import Graph, URIRef, Namespace
from rdflib.namespace import SH, RDF

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
ontologyPath = os.path.join(baseDir2,'ConcreteOntology')
metadataPath = os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv')
graphPath = os.path.join(baseDir0,'E-modul-processed-data/EM_Graph.ttl')
#importedOntologiesPath = os.path.join(baseDir2,'GraphCreation/Ontologies')
processedDataPath = os.path.join(baseDir0,'E-modul-processed-data')

SCHEMA = Namespace('http://schema.org/')

"""
Given a path to a shacl shape and a path to an rdf file, this function tests the rdf data against the specified shacl shapes.
The result is an rdflib graph containing the validation report, if it is empty the validation was successful.
"""
def test_file(shapepath: str, filepath: str) -> Graph:
    with open(shapepath, 'r') as f:
        shacl_shape = f.read()
    with open(filepath, 'r') as f:
        rdf_data = f.read()

    # post both files to the validation endpoint of the docker app
    post_data = {'rdf': rdf_data, 'shaclShape': shacl_shape}
    r = s.post('http://localhost:8080/rdfconv/validate', files=post_data)
    result_graph = Graph()

    # only parse graph if it contains any violations
    if not r.text.startswith('VALID'):
        res = r.text[6:] #remove the 'fail: ' prefix from the result
        result_graph.parse(data=res, format='ttl')

        # also add nodes from data and shacl shapes to graph to be able to search backwards for the violated shapes
        shacl_graph = Graph()
        data_graph = Graph()
        shacl_graph.parse(data=shacl_shape, format='ttl')
        data_graph.parse(data=rdf_data, format='ttl')
        result_graph += shacl_graph
        result_graph += data_graph
    
    return result_graph

"""
Returns true if the given shape is violated in the report.
"""
def violates_shape(validation_report: Graph, shape: URIRef):

    # get the class that is targeted by the specified shape
    target_class = validation_report.value(shape, SH.targetClass, None, any=False)
    if target_class is None:
        raise ValueError(f'The shapes graph does not contain a {shape} shape.')

    
    # get all classes that have been violated
    # check if any of the violated classes is the class that is targeted by the specified shape
    for o in validation_report.objects(None, SH.focusNode):
        if target_class in validation_report.objects(o, RDF.type, None):
            return True

    # no violated class is targeted by the specified shape, thus the shape is not violated
    return False


"""
Clone and run the RDFConverter microservice through docker. This allows access to the schacl validation api via POST request.
"""
Repo.clone_from('https://github.com/Mat-O-Lab/RDFConverter.git', os.path.join(baseDir1, 'RDFConverter'))
subprocess.run(['docker-compose', '-f', os.path.join(baseDir1, 'RDFConverter/docker-compose.yml'), 'up', '-d', 'app'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

# read shacl shape
shacl_path = os.path.join(baseDir1, 'shape.ttl')
sh_graph = Graph()

# sleep shortly to avoid access of the validation service before it is ready
sleep(3)
s = requests.session()

# test the EM_Graph.ttl graph that was created in the steps before against the shacl shape
rdf_file = os.path.join(processedDataPath, 'EM_Graph.ttl')
g = test_file(shacl_path, rdf_file)

# close session object and shutdown docker container
s.close()
subprocess.run(['docker-compose', '-f', os.path.join(baseDir1, 'RDFConverter/docker-compose.yml'), 'down'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


# Tests specific assertions about the shacl validation report.
# assert that certain violations did not occur:
assert not violates_shape(g, SCHEMA.SpecimenDiameterShape)
assert not violates_shape(g, SCHEMA.SpecimenShape)

# assert that certain violations occurred:
assert violates_shape(g, SCHEMA.InformationBearingEntityShape)
