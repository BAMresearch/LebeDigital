## News/Info
- The KG's now use the concrete ontology based on the new pmd core ontology (https://github.com/materialdigital/core-ontology/blob/main/pmdco_core.ttl) and the QUDT Ontology for the units (https://github.com/qudt/qudt-public-repo/blob/main/vocab/unit/VOCAB_QUDT-UNITS-ALL-v2.1.ttl)
- Moved all Base-Ontologies in the BaseOntology folder
- The Concrete Ontology is still work in progress

## Knowledge Graphs Desciption

- CPTO: The Concrete Production and Testing Ontology (CPTO) is created in the project LeBeDigital.
        It consists of the: Heat of Hydration Ontology, Mixture Design Ontology, Mixing Process Ontology,
        Hardening Ontology, Compressive Strength Ontology and the Secant Modulus of Elasticity Ontology.
- Heat of Hydration Ontology: Describes the process for determining the Heat of Hydration.
- Mixture Design Ontology: Describes the ingridients and mixing ratios of the concrete used. Also contains the mixing date and location.
- Mixing Process Ontology: Describes the concrete-mixing process e.g. which mixer was used.
- Hardening Ontology: Describes the conditions while the concrete is hardening.
- Compressive Strength Ontology: Describes the compressive strength test of each sample.
- Secant Modulus of Elasticity Ontology: Describes testing the secant elastic modulus of each sample.

## Documentation:

### Exporting the Knowledge Graphs to ttl using Ontopanel

- Open the diagrams.net file using the following link: https://yuechenbam.github.io/src/main/webapp/index.html (Accept the Prompt for the Plugins, if it doesn't 
show when opening go to Extra -> Plugins -> Add... then add importOnto and ontopanel. After that reload the page)
- (Change the Ontology etc.)
- To Export: go to Extras -> Ontopanel-Converter
- Choose the format you need then click on Convert
- Go to Show Error to check for errors (Ignore the No data for mapping error)
- Click on the red downwards arrow to download the file

### The Concrete Ontology (Base-Ontology for all Knowledge Graphs)

- The Concrete Ontology is found in the folder "BaseOntology"
- The Ontology includes all taxonomie described in the excel file "Taxo.xlsx"
- The tool tax2ont was used for developing the Ontology and could be useful in the future, further information is in the readme there.
- The Concrete Ontology is using the current pmd core ontology and the QUDT Ontology for the Units, a "dictionary" for the mapping of the units can be found on (https://www.qudt.org/doc/DOC_VOCAB-UNITS.html)