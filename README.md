## Introction - Model calibration
Model calibration (or model updating) is a task most scientists are facing when building parameterized models and then
 updating the parameters based on some experimental data in order to generalize the model and make accurate predictions.
For complex models and complex data sets (different experimental tests with different data structures, sensor
types and different models) this task is often tedious, difficult to reproduce and often error prone due to complex
challenges related to data processing, forward model development and inference. The aim of this roject is to make
the process more transparent, easier to setup and work with and more transparent when different people are jointly
performing this task.
Further information can be found in the [documentation](https://modelcalibration.readthedocs.io/en/latest/?).

## Installation of the conda environment
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


