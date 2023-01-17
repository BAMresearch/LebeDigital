import os
from pathlib import Path

from lebedigital.shacl.validation import SCHEMA, read_graph_from_file, test_graph as graph_test, violates_shape

def test_graph_against_shacl_shape():
    """
    Testing a data graph against a shapes graph and checking that the appropriate shapes fail.
    """
    shacl_directory = os.path.dirname(Path(__file__))
    shapes_location = Path(shacl_directory, 'shape.ttl')
    data_location = Path(shacl_directory, 'graph.ttl')

    g = read_graph_from_file(data_location)
    s = read_graph_from_file(shapes_location)

    res = graph_test(g, s)

    assert not violates_shape(res, SCHEMA.SpecimenDiameterShape)
    assert not violates_shape(res, SCHEMA.SpecimenShape)
    assert violates_shape(res, SCHEMA.InformationBearingEntityShape)

