# Installation conda environment
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

#General Remarks
The example solves a problem with a virtual model in specified in Correlation.py. The functions is 
defined by
```
f(x) = a+b*x+c*x^2.
```
Measurements are performed at a given number of sensor locations
```
x_function_sensors = np.linspace(0, 1, options['n_function_sensors'])
x_derivative_sensors = np.linspace(0, 1, options['n_derivative_sensors']).
```
Thus you can create function sensors and sensors that measure the derivative of the function
at different locations in the interval [0,1].

All options are specified in a yaml file in the directory. If you want to add another setup, just cp a yaml file, rename it.  
The computation is initiated via
```
doit
```
This is a [toolbox in python](https://pydoit.org/) (similar to make, or other workflow environments that allows define dependencies between your tasks and only update the dependent task if you have changed something in your tool chain.)
More information can e.g. be found [here](https://bamresearch.github.io/Reproducible-Science/).

The definition of the variables in a yaml file are e.g.
```
n_function_sensors: 50
n_derivative_sensors: 0
n_configurations: 3

correlation_length_function: 0.3
correlation_length_derivative: 0.3

virtual_experiment_correlation_length_function: 0.3
virtual_experiment_correlation_length_derivative: 0.3

virtual_experiment_noise_std_function: 0.3
virtual_experiment_noise_std_derivative: 0.3

virtual_experiment_b: 3
virtual_experiment_c: 2.5
```
The number of configurations specifies different configurations represent by the parameter a
```
a = np.linspace(0,4,options['n_configurations']
```
so a is a value in the interval [0,4]. Measurements are then given for each a, thus the total number of 
samples is n_configurations*(n_function_sensors+n_derivative_sensors)

The options correlation_length_function and correlation_length_derivative specify the correlation 
length for the noise (for the function and derivative sensors separately). 

The options for the (virtual) experiment that is used to generate the data are given at the end. Note that c 
is a quadratic term that is not part of the function to be estimated.
```
f(x) = a+b*x
f_virt(x) = a+b*x+c*x^2
```
This term is not present in the model to identified and represents the model bias (the part of the data that is not explainable) by the model.
The noise in the virtual experiment is also assumed to be correlated. If you want to remove the correlation, simply set the correlation
length to zero. 

The output for each yaml file is a png file with the virtual data points, the exact virtual function as well as the approximation
with the predictive posterior including the uncertainty related to the identified noise term. The parameters of the posterior are also 
given in the command line.