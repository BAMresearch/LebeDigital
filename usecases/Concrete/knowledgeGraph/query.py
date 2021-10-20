import rdflib
from SPARQLWrapper import SPARQLWrapper
from pathlib import Path
import os
import sys

baseDir1 = Path(__file__).resolve().parents[1]
baseDir0 = Path(__file__).resolve().parents[0]
triplePath = os.path.join(baseDir0,'E-modul-processed-data/EM_Graph.ttl')
print(baseDir0)
if sys.platform == 'win32':
    prefixPath = 'file:///' + os.path.join(baseDir0,'E-modul-processed-data').replace('\\','/') + '/'
else:
    prefixPath = 'file://' + os.path.join(baseDir0,'E-modul-processed-data') + '/'

graph = rdflib.Graph()
graph.parse(triplePath, format='n3')

print('---------------------------------------------------------------------')
print(prefixPath)
print('sample queries')

q = f"""
    prefix bwmd: <{prefixPath}https%3A//www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#>
    prefix mseo: <{prefixPath}https%3A//purl.matolab.org/mseo/mid/>
    prefix cco: <{prefixPath}http%3A//www.ontologyrepository.com/CommonCoreOntologies/>
    prefix obo: <{prefixPath}http%3A//purl.obolibrary.org/obo/>
    prefix con: <{prefixPath}https%3A//github.com/BAMresearch/ModelCalibration/blob/Datasets/usecases/Concrete/ConcreteOntology/Concrete_Ontology_MSEO.owl#>
    select ?rawdatapath
    where {{
        {{
            select ?rawdatafile
            where {{
                {{
                    select ?outputfile
                    where {{
                        mseo:E-modul_experiment_BA_Los_M_V-4 cco:has_output ?outputfile
                    }}
                }}
                ?outputfile
                obo:RO_0010001
                ?rawdatafile
            }}
        }}
        ?rawdatafile
        cco:has_URI_value
        ?rawdatapath
    }}
        
"""
q1 = """
    select ?s ?p ?o
    where {
        ?s ?p ?o
    }
"""

results = graph.query(q)
for result in results:
    if sys.platform == 'win32':
        print(f"{result}".encode("utf-8"))
    else:
        print(result)

