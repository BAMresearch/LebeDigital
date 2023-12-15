## Info
- The Knowledge Graphs use the Concrete Ontology (https://w3id.org/cpto) based on the new Pmd core ontology (https://github.com/materialdigital/core-ontology/blob/main/pmdco_core.ttl) and the QUDT Ontology for the units (https://github.com/qudt/qudt-public-repo/blob/main/vocab/unit/VOCAB_QUDT-UNITS-ALL-v2.1.ttl)
- The Base Ontology (CPTO - Concrete Production and Testing Ontology) is located in the base Ontology folder and can be accessed via https://w3id.org/cpto

## Knowledge Graphs Description

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

### The w3id link (https://w3id.org/cpto)

- Persistent URI Service: w3id.org offers a persistent and stable Uniform Resource Identifier (URI) service, primarily for the Semantic Web community
- Avoids Broken Links: Its main goal is to prevent broken links and ensure long-term accessibility of web resources, such as ontologies or vocabularies used in Semantic Web and Linked Data projects
- Permanent Redirects: w3id.org creates permanent HTTP URIs that redirect to actual content locations. This redirection can be updated if the content's location changes. Currently it redirects to the raw page on github (https://raw.githubusercontent.com/BAMresearch/LebeDigital/main/lebedigital/ConcreteOntology/BaseOntology/ConcreteOntology.owl)
- GitHub Repository: Changes to w3id.org URIs are managed through their GitHub repository (https://github.com/perma-id/w3id.org). Users submit change requests via this repository
- Pull Request Process: To modify a link, you need to submit a pull request with the proposed changes to the w3id.org GitHub repository
- Required Information: The request should include the new target URL and possibly a justification for the change
- Community Review: The change is reviewed by the community maintainers. If approved, the URI will redirect to the new target URL

### The Concrete Ontology (Base-Ontology for all Knowledge Graphs)

- The Concrete Ontology is found in the folder "BaseOntology"
- The Ontology includes all info described in the excel file "Taxonomie.xlsx"
- The Concrete Ontology is using the current PMD-Core ontology and the QUDT Ontology for the Units, a "dictionary" for the mapping of the units can be found on (https://www.qudt.org/doc/DOC_VOCAB-UNITS.html)
