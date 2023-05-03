[![Continuous integration](https://github.com/BAMresearch/LebeDigital/actions/workflows/lebedigital.yml/badge.svg)](https://github.com/BAMresearch/LebeDigital/actions)
[![coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/eriktamsen/c10a5b6d0714b1fe2344eb60918e92f8/raw/lebedigital_main_coverage.json)](https://en.wikipedia.org/wiki/Code_coverage)


## Overview
This repository contains an automatic workflow to extract relevant 
metadata from concrete tests (e.g. Young's modulus, compression tests).  

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
