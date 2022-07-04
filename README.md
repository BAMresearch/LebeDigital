## Overview
This repository contains an automatic workflow to extract relevant 
metadata from concrete tests (e.g. Young's modulus, compression tests).  

## Installation of the conda environment
Clone the repository and make sure to include the submodules
```
git clone https://github.com/BAMresearch/ModelCalibration.git
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