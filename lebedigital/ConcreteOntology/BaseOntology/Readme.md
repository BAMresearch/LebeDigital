# This is a test for creation an ontology using owlready
[Owlready](https://owlready2.readthedocs.io) is a tool for generating ontologies and also knowledge bases. 

## Installation
Owlready2 is a python library that can be installed using pip or conda and could be just added to the environment by either 
adding it to the environment.yml, or in conda using the command:
```bash
conda install -c conda-forge owlready2
``` 
## Using external ontologies (pmdco, qudt)
Owlready only works with a limited number of formats (RDF/XML, OWL/XML, NTriples), thus some ontologies have to be
converted. The [PMDcore ontology](https://github.com/materialdigital/core-ontology) can be accessed via its w3id
[https://w3id.org/pmd/co](https://w3id.org/pmd/co), which directly provides a 
[link](https://materialdigital.github.io/core-ontology/ontology.rdf) to the ontology in RDF/XML format. A local copy is 
stored as well (co.rdf), where some bugs are corrected (an [issue](https://github.com/materialdigital/core-ontology/issues/67)
 is filed to correct the bug. Other ontologies (e.g. [QUDT](http://qudt.org/2.1/schema/qudt)) can be downloaded
and then converted to RDF/XML using e.g. [owl-cli](https://atextor.de/owl-cli/). 
Follow the installation instructions, download the qudt and then use the command:
```bash
 ./owl-x86_64-linux-snapshot  write -o rdfxml SCHEMA_QUDT-v2.1.ttl >& qudt.rdf
```
This can then also be used to convert the created CPTO from the RDF/XML format to turtle or other formats.
```bash
./owl-x86_64-linux-snapshot write -i rdfxml -o turtle CPTO.rdf &> CPTO.ttl
```