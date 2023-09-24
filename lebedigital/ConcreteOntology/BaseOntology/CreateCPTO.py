from owlready2 import *

cpto = get_ontology("https://w3id.org/cpto/co#")

pmdcore = get_ontology(
    "file:///home/junger/github/LebeDigital/lebedigital/ConcreteOntology/BaseOntology/co.rdf"
).load()

# converted using https://atextor.de/owl-cli/
# ./owl-x86_64-linux-snapshot  write -o rdfxml SCHEMA_QUDT-v2.1.ttl >& SCHEMA_QUDT-v2.1.rdf
qudt = get_ontology("file:///home/junger/github/LebeDigital/lebedigital/ConcreteOntology/BaseOntology/qudt.rdf").load()

#cpto.imported_ontologies.append(pmdcore)
#cpto.imported_ontologies.append(qudt)

with cpto:
    class ModulusOfElasticity(pmdcore.ProcessingNode):
        pass
    class Specimen(pmdcore.TestPieceInfo):
        pass


ModulusOfElasticity.comment = ""
ModulusOfElasticity.label = ['"Bestimmung des Elastizitaetsmoduls (Sekantenmodul)"@de', '"Determination of Secant Modulus of Elasticity"@en']

Specimen.comment = ""
Specimen.label = ['"Probenalter"@de', '"specimen age"@en']


myE = ModulusOfElasticity("myURI_E")
mySpecimen = Specimen("myURI_Specimen")
pmdcore.input.python_name = "input"
myE.input = [mySpecimen]

cpto.save("CPTO.rdf", format="rdfxml")





