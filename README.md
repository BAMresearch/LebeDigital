[![Continuous integration](https://github.com/BAMresearch/LebeDigital/actions/workflows/lebedigital.yml/badge.svg)](https://github.com/BAMresearch/LebeDigital/actions)
[![DOI](https://zenodo.org/badge/311899936.svg)](https://zenodo.org/badge/latestdoi/311899936)


## Overview
The term “LeBe” from LeBeDigital stands in German for “Lebenszyklus von Beton” and in English for “Concrete life cycle”. LeBeDigital is about developing an ontology for the concrete production process chain. The Project LeBeDigital is sponsored by BMBF and runs under the supervision of the Plattform MaterialDigital. The LeBeDigital-team consist of 3 cooperations partners: The Federal Institute for Materials Research and Testing (BAM), the Karlsruhe Institut of Technologie (KIT) and the Technical University of Munich (TUM).

Concrete, the main subject of the project, is a very complex material considering the diversity of the base materials and component designs, the high complexity of recipes and the manufacturing process, the diversity in molding designs as well as the time dependency of all properties.
The expected outcome of the project is to develop a material database, where concrete-specific characteristic values and models are structurally integrated. The main objectives are:
1. Ontology-based collection of data and metadata in TripleStores
2. Linking existing, experience-based, experimental and simulated data
3. Combined material optimisation and geometry optimisation of precast elements based on simulation workflows
4. Model calibration

## Installation of the conda environment
Clone the repository and make sure to include the submodules
```
git clone https://github.com/BAMresearch/LebeDigital.git
```

Create environment within the main project folder
```
conda env create --prefix ./conda-env -f environment.yml 
```

Update environment within main project folder (if environment.yml was changed)
```
conda env update --prefix ./conda-env -f environment.yml --prune
```

Activate environment
```
conda activate ./conda-env
```

## Structure
All the implemented functions that perform a certain action should be 
implemented in a local modul called lebedigital. This is installed thus 
allowing to import these classes and functions in other files, e.g. the test 
files or the usecases. For the usecases, we intend to have two, one is a 
minimum working example, and another one relates to the final demonstrator. 
