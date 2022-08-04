from pyshacl import validate
from rdflib import Graph, URIRef, Namespace
from rdflib.util import guess_format
from rdflib.namespace import SH, RDF

SCHEMA = Namespace('http://schema.org/')


def test_graph(rdf_graph: Graph, shapes_graph: Graph) -> Graph:
    """
    Tests an RDF graph against a SHACL shapes graph.

    Parameters
    ----------
    rdf_graph
        An rdflib Graph object containing the triples to test against.
    shapes_graph
        An rdflib Graph object containing the shapes to test.

    Returns
    -------
    result_graph
        An rdflib Graph object containing the SHACL validation report (which is empty if no SHACl shapes were violated).
    """
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

def violates_shape(validation_report: Graph, shape: URIRef) -> bool:
    """
    Returns true if the given shape is violated in the report.

    Parameters
    ----------
    validation_report
        An rdflib Graph object containing a validation report from the test_graph function.
    shape
        A URIRef object containing the URI of a shape.
    
    Returns
    -------
        True, if the specified shape appears as violated in the validation report, False otherwise.
    """
    # get the class that is targeted by the specified shape
    target_class = validation_report.value(shape, SH.targetClass, None, any=False)
    if target_class is None:
        raise ValueError(f'The shapes graph does not contain a {shape} shape.')

    
    # get all classes that have been violated
    # check if any of the violated classes is the class that is targeted by the specified shape
    # return any((True for o in validation_report.objects(None, SH.focusNode) if target_class in validation_report.objects(o, RDF.type)))
    for o in validation_report.objects(None, SH.focusNode):
        if target_class in validation_report.objects(o, RDF.type):
            return True

    # no violated class is targeted by the specified shape, thus the shape is not violated
    return False


def read_graph_from_file(filepath: str) -> Graph:
    """
    Reads a file containing an RDF graph into an rdflib Graph object.

    Parameters
    ----------
    filepath
        The path to the file containing the graph.

    Returns
    -------
    graph
        The rdflib Graph object containing the triples from the file.
    """
    with open(filepath, 'r') as f:
        graph = Graph()
        graph.parse(file=f, format=guess_format(filepath))
    return graph