from pyshacl import validate
from rdflib import Graph, URIRef, Namespace
from rdflib.util import guess_format
from rdflib.namespace import SH, RDF

"""
baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
ontologyPath = os.path.join(baseDir2,'ConcreteOntology')
metadataPath = os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv')
graphPath = os.path.join(baseDir0,'E-modul-processed-data/EM_Graph.ttl')
processedDataPath = os.path.join(baseDir0,'E-modul-processed-data')
"""

SCHEMA = Namespace('http://schema.org/')

"""
Given a path to a shacl shape and a path to an rdf file, this function tests the rdf data against the specified shacl shapes.
The result is an rdflib graph containing the validation report, if it is empty the validation was successful.
"""
def test_graph(rdf_graph: Graph, shapes_graph: Graph) -> Graph:

    conforms, result_graph, _ = validate(
            rdf_graph,
            shapes_graph,
            ont_graph=None,  # can use a Web URL for a graph containing extra ontological information
            inference='none',
            abort_on_first=False,
            allow_infos=False,
            allow_warnings=False,
            meta_shacl=False,
            advanced=False,
            js=False,
            debug=False)

    # only add other graphs if any violations occurred
    if not conforms:
        # also add nodes from data and shacl shapes to graph to be able to search backwards for the violated shapes
        result_graph += shapes_graph
        result_graph += rdf_graph
    
    return result_graph

"""
Returns true if the given shape is violated in the report.
"""
def violates_shape(validation_report: Graph, shape: URIRef) -> bool:

    # get the class that is targeted by the specified shape
    target_class = validation_report.value(shape, SH.targetClass, None, any=False)
    if target_class is None:
        raise ValueError(f'The shapes graph does not contain a {shape} shape.')

    
    # get all classes that have been violated
    # check if any of the violated classes is the class that is targeted by the specified shape
    for o in validation_report.objects(None, SH.focusNode):
        if target_class in validation_report.objects(o, RDF.type):
            return True

    # no violated class is targeted by the specified shape, thus the shape is not violated
    return False

"""
Reads a graph from a file into a Graph object.
"""
def read_graph_from_file(filepath: str) -> Graph:
    with open(filepath, 'r') as f:
        graph = Graph()
        graph.parse(file=f, format=guess_format(filepath))
    return graph


# assert that certain violations occurred / did not occur:
# assert violates_shape(g, SCHEMA.InformationBearingEntityShape)
# assert not violates_shape(g, SCHEMA.InformationBearingEntityShape)

