from owlready2 import *
from pathlib import Path

cpto = get_ontology("https://w3id.org/cpto/co#")

pmdcore_path = "file://" + Path(Path(__file__).parent.resolve() / "co.rdf").as_posix()
pmdcore = get_ontology(pmdcore_path).load()

qudt_path = "file://" + Path(Path(__file__).parent.resolve() / "qudt.rdf").as_posix()
qudt = get_ontology(qudt_path).load()

#cpto.imported_ontologies.append(pmdcore)
#cpto.imported_ontologies.append(qudt)

with cpto:
    class ModulusOfElasticity(pmdcore.ProcessingNode):
        pass
    ModulusOfElasticity.comment = ""
    ModulusOfElasticity.label = [locstr("Bestimmung des Elastizitaetsmoduls (Sekantenmodul)", lang="de"),
                                 locstr("Determination of Secant Modulus of Elasticity", lang="en")]
    class Specimen(pmdcore.TestPieceInfo):
        pass
    Specimen.comment = ""
    Specimen.label = ['"Probenalter"@de', '"specimen age"@en']





#myE = ModulusOfElasticity("myURI_E")
#mySpecimen = Specimen("myURI_Specimen")
#pmdcore.input.python_name = "input"
#myE.input = [mySpecimen]

cpto.save("CPTO.rdf", format="rdfxml")





