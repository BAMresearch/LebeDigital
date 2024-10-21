<h1 align="center">LeBeDigital</h1>
<h3 align="center">Life Cycle of Concrete</h3>


<p align="center">
  <a href="#dart-project-objective">Objective</a> &#xa0; | &#xa0; 
  <a href="#sparkles-key-features">Features</a> &#xa0; | &#xa0;
  <a href="#rocket-technologies">Technologies</a> &#xa0; | &#xa0;
  <a href="#white_check_mark-installation">Installation</a> &#xa0; | &#xa0;
  <a href="#memo-license">License</a> 
</p>

<br>

## :dart: Project Objective ##

The aim of the joint project LeBeDigital, a part of the initiative Plattform Material Digital, is to develop a structure for concrete material data following the FAIR principle and being flexible enough to fulfill requirements of future developments in concrete research and production. The expected outcome is a material database, where Concrete-specific characteristic values and models are structurally integrated.

## :sparkles: Key Features ##

:heavy_check_mark: Raw Data Upload\
:heavy_check_mark: Metadata Extraction\
:heavy_check_mark: Knowledge Graph Mapping\
:heavy_check_mark: Data Retrieval

## :rocket: Technologies ##

- [Python](https://www.python.org)
- [Flask](https://flask.palletsprojects.com/)
- [PMD Core Ontology](https://github.com/materialdigital/core-ontology)
- [SPARQL](https://www.w3.org/TR/sparql11-query/)
- [Ontodocker](https://materialdigital.github.io/pmd-server/pages/services/onto-docker/)


## :white_check_mark: Installation ##

Before starting, you need to have [Git](https://git-scm.com) and [Python](https://www.python.org) installed.


```bash
# Clone this project
$ git clone https://github.com/BAMresearch/LebeDigital.git

# Change branch
$ git checkout workflowTest 

# Activate virtual environment
# For Linux:
$ source venv/bin/activate 

# For Windows:
$ .\venv\Scripts\activate 
    
# Access
$ cd web

# Install dependencies
$ pip install -r requirements.txt

# Set up configuration

$ cp config.template.json config.json 
# Open `config.json` and replace placeholders with actual credentials.

# Run the project
$ python server.py

```

## :memo: License ##

This project is under license from MIT. For more details, see the [LICENSE](LICENSE) file.


<a href="#top">Back to top</a>
