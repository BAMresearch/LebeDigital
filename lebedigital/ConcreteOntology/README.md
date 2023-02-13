## Info
- For each Ontology there is a drawio-file and a ttl file
- The Mapping Example is moved to usecases/MinimumWorkingExample/Example_Mapping
- Moved all Base-Ontologies in the BaseOntology folder

## Ontologies:

- CPTO: The Concrete Production and Testing Ontology (CPTO) is created in the project LeBeDigital.
        It consists of the: Heat of Hydration Ontology, Mixture Design Ontology, Mixing Process Ontology,
        Hardening Ontology, Compressive Strength Ontology and the Secant Modulus of Elasticity Ontology.
- Heat of Hydration Ontology: Describes the process for determining the Heat of Hydration.
- Mixture Design Ontology: Describes the ingridients and mixing ratios of the concrete used. Also contains the mixing date and location.
- Mixing Process Ontology: Describes the concrete-mixing process e.g. which mixer was used.
- Hardening Ontology: Describes the conditions while the concrete is hardening.
- Compressive Strength Ontology: Describes the compressive strength test of each sample.
- Secant Modulus of Elasticity Ontology: Describes testing the secant elastic modulus of each sample.

## Exporting the Ontology to ttl using Ontopanel

- Open the diagrams.net file using the following link: https://yuechenbam.github.io/src/main/webapp/index.html (Accept the Prompt for the Plugins, if it doesn't 
show when opening go to Extra -> Plugins -> Add... then add importOnto and ontopanel. After that reload the page)
- Change the Ontology etc.
- To Export: go to Extras -> Ontopanel-Converter
- Choose the format you need then click on Convert
- Go to Show Error to check for errors (Ignore the No data for mapping error)
- Click on the Red downwards arrow to Download the file
