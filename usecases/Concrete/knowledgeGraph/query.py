import rdflib
from SPARQLWrapper import SPARQLWrapper
from pathlib import Path
import os

baseDir1 = Path(__file__).resolve().parents[1]
baseDir0 = Path(__file__).resolve().parents[0]
triplePath = os.path.join(baseDir0,'E-modul-processed-data/EM_Graph.ttl')

graph = rdflib.Graph()
graph.parse(triplePath, format='n3')

print('---------------------------------------------------------------------')
print('sample queries')

q = """
    prefix a: <http://www.ontologyrepository.com/CommonCoreOntologies/>
    select ?s ?p ?o
    where {
        ?s ?p "Kh"
    }
    limit 10
"""

results = graph.query(q)
for result in results:
    print(result)