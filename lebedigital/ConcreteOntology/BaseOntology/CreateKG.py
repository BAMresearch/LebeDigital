from owlready2 import *
from pathlib import Path

export = get_ontology("LebeDigitalData")

pmdcore_path = "file://" + Path(Path(__file__).parent.resolve() / "co.rdf").as_posix()
pmdcore = get_ontology(pmdcore_path).load()

#qudt_path = "file://" + Path(Path(__file__).parent.resolve() / "unit.ntriple").as_posix()
#qudt = get_ontology(qudt_path).load()

cpto_path = "file://" + Path(Path(__file__).parent.resolve() / "cpto.rdf").as_posix()
cpto = get_ontology(cpto_path).load()

export.imported_ontologies.append(cpto)

# define python names for properties from PMD core
pmdcore.input.python_name = "input"
pmdcore.characteristic.python_name = "characteristic"
pmdcore.value.python_name = "value"
pmdcore.unit.python_name = "unit"

mySpecimen = cpto.Specimen("myURI_Specimen", namespace=export)
my_E = cpto.ModulusOfElasticity("myURI_E", namespace=export)
my_diameter = pmdcore.Diameter("my_diameter", namespace=export)
my_diameter.value = [10]
#my_diameter.unit = [qudt.MilliM]

my_E.characteristic = [my_diameter]
my_E.input = [mySpecimen]

export.save("myKG.rdf", format="rdfxml")
