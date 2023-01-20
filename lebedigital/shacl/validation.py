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
    _, result_graph, _ = validate(
            rdf_graph,
            shacl_graph=shapes_graph,
            ont_graph=None,  # can use a Web URL for a graph containing extra ontological information
            inference='none',
            abort_on_first=False,
            allow_infos=False,
            allow_warnings=False,
            meta_shacl=False,
            advanced=False,
            js=False,
            debug=False)

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
    if (shape, RDF.type, None) not in validation_report:
        raise ValueError(f'The shacl shape graph does not contain a {shape} shape.')

    # get the properties that are used by the specified shape
    properties = validation_report.triples((shape, SH.property, None))

    # check if any property of the shape appears in the validation report (then it was violated)
    for s, p, o in properties:
        if (None, SH.sourceShape, o) in validation_report:
            return True;

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

def violates_shapes_list(res: Graph, shapes_list: list[URIRef]) -> bool:
    """
    Reads a data and a shapes graph from files and tests a list of shapes from the shapes graph against the data graph.

    Parameters
    ----------
    graph_path
        The data graph.
    shapes_graph_path
        The shapes graph.
    shapes_list
        The list of shapes from the shapes graph that should be tested on the data graph.

    Returns
    -------
        True, if the graph violates any of the shapes from the shapes list.
    """
    
    return any(violates_shape(res, shape) for shape in shapes_list)