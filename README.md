# ModelCalibration
Development of general interfaces to make model calibration tasks transparent

# Installation of the conda environment
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

Install submodule Bayes
```
git pull --recurse-submodules
pip install -e BayesianInference
```


