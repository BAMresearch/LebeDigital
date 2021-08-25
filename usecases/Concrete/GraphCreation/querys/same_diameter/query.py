from rdflib import Graph, URIRef, Literal
from rdflib.namespace import RDFS
from SPARQLWrapper import SPARQLWrapper, JSON, POST, BASIC, DIGEST
import yaml

data = {}
query = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX cc: <http://creativecommons.org/ns#>
prefix mid: <https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#>
prefix WCTmid: <https://mobi.com/ontologies/6/2021/WCTmid#>
prefix WCT: <https://mobi.com/ontologies/6/2021/WoodCompressionTest#>
prefix cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
prefix obo: <http://purl.obolibrary.org/obo/>
prefix cst: <https://mobi.com/ontologies/7/2021/ConcreteStressTestOntologie#>

SELECT ?s ?height ?heightvalue ?heightICE ?heightIB ?heightUnit ?diameter ?diametervalue ?diameterICE  ?diameterIB ?diameterUnit ?weight ?weightvalue ?weightICE ?weightIB ?weightUnit ?rezeptICE ?rezeptIB ?rezept_uri ?loaddisplacementdata
WHERE {
  ?s rdf:type cst:ConcreteCompressionTest.
  ?Specimen cco:is_affected_by ?s.
  ?output cco:is_output_of ?s.

  ?output obo:RO_0010001 ?outputIB.
  ?outputIB cco:has_URI_value ?loaddisplacementdata.

  ?Specimen obo:RO_0000086 ?weight.
  ?weight rdf:type cco:Weight.
  ?weight cco:is_measured_by ?weightICE.
  ?weightICE obo:RO_0010001 ?weightIB.
  ?weightIB cco:has_decimal_value ?weightvalue.
  ?weightIB cco:uses_measurement_unit ?weightUnit.


  ?Specimen obo:RO_0000086 ?diameter.
  ?diameter rdf:type cco:Diameter.
  ?diameter cco:is_measured_by ?diameterICE.
  ?diameterICE obo:RO_0010001 ?diameterIB.
  ?diameterIB cco:has_decimal_value ?diametervalue.
  ?diameterIB cco:uses_measurement_unit ?diameterUnit.

  ?Specimen obo:RO_0000086 ?height.
  ?height rdf:type cco:Height.
  ?height cco:is_measured_by ?heightICE.
  ?heightICE obo:RO_0010001 ?heightIB.
  ?heightIB cco:has_decimal_value ?heightvalue.
  ?heightIB cco:uses_measurement_unit ?heightUnit.


}
"""
###############################################################################

print('for querying an dataset stored at "https://matolab.bam.de", an connection must be established')
print('make sure you are already connected via vpn to "cvpn.bam.de"')
print('establishing connection to https://matolab.bam.de')
username = input("enter your username: ")
password = input("enter your password: ")
graph_name = input("enter the name of the dataset you want to query (this name must match perfektly to the name of an existing one): ")


sparql = SPARQLWrapper("https://matolab.bam.de/graph/{}/sparql".format(graph_name))
sparql.setHTTPAuth(BASIC)
sparql.setCredentials(username, password)
sparql.setMethod(POST)

sparql.setQuery(query)

results = sparql.query()

sparql.setReturnFormat(JSON)
results = sparql.query().convert()
i = 0
Diameter_list = [100,200,300]
for result in results["results"]["bindings"]:

    for diameter in Diameter_list:
        if diameter -10 <= float(result['diametervalue']['value']) <= diameter+10:

            if '{}'.format(diameter) not in data:
                data['{}'.format(diameter)] = {}



            data['{}'.format(diameter)]['Experiment{}'.format(i)] = {}
            data['{}'.format(diameter)]['Experiment{}'.format(i)]['height'] = {'value' : result['heightvalue']['value'] , 'unit' : result['heightUnit']['value'].split('/')[-1]}
            data['{}'.format(diameter)]['Experiment{}'.format(i)]['diameter'] = {'value' : result['diametervalue']['value'] , 'unit' : result['diameterUnit']['value'].split('/')[-1]}
            data['{}'.format(diameter)]['Experiment{}'.format(i)]['weight'] = {'value' : result['weightvalue']['value'] , 'unit' : result['weightUnit']['value'].split('/')[-1]}
            data['{}'.format(diameter)]['Experiment{}'.format(i)]['LoadDisplacementData'] = result['loaddisplacementdata']['value']

    i+=1
print(data)
with open('Experiments_with_same_diameters.yml', 'w') as f:
    yaml.dump(data, f)
