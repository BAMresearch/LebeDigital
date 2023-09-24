from owlready2 import *

cpto = get_ontology(
    "file:///home/junger/github/LebeDigital/lebedigital/ConcreteOntology/BaseOntology/CPTO.rdf"
).load()

pmdcore = get_ontology(
    "file:///home/junger/github/LebeDigital/lebedigital/ConcreteOntology/BaseOntology/co.rdf"
).load()


myE = cpto.ModulusOfElasticity("myURI_E")
mySpecimen = cpto.Specimen("myURI_Specimen")
pmdcore.input.python_name = "input"
myE.input = [mySpecimen]

cpto.save("myKG.rdf", format="rdfxml")
