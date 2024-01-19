"""
Get version of an ontology
"""
import argparse
import rdflib

parser = argparse.ArgumentParser(prog='get_ontoversion', description='Get version of an ontology')
parser.add_argument('ontology_file', help='Location of the ontology file')
args = parser.parse_args()

g = rdflib.Graph()
g.parse(location=args.ontology_file)

# Get first URI and version of the ontology (assuming there is only one)
QUERY_OBJECT="""
PREFIX owl: <http://www.w3.org/2002/07/owl#>
SELECT ?s ?o WHERE {
    ?s a owl:Ontology ;
       owl:versionInfo ?o .
}
LIMIT 1
"""
for row in g.query(QUERY_OBJECT):
    print(row.o)
    break
