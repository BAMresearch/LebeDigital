import numpy as np
import yaml


def write_metadata_to_yaml(virtual_experimental_metadata_file,
                           num_a, b, c,
                           num_function_sensors, sigma_noise_function,
                           num_derivative_sensors, sigma_noise_derivative,
                           center=False, seed=42): 	
    """Create meta data of virtual experiment
    Args:
        virtual_experimental_metadata_file: filename for metadata to be generated
        num_a: number of constant offset a computed by 1+2*i with i in range(num_a)
        b: linear coefficient of the model
        c: quadratic coefficient of the model
        num_function_sensors(int): number of derivative sensors in the interval [0,1]
        sigma_noise_function: std deviation of the function noise
        num_derivative_sensors(int): number of derivative sensors in the interval [0,1]
        sigma_noise_derivative: std deviation of the derivative noise
        center: if False, include start and end point, if True, divide in equal intervals and use midpoint
        seed: seed for the random number generator
    """

    if center is False:
        x_function = np.linspace(70, 90, num_function_sensors)
        x_derivative = np.linspace(0, 1, num_derivative_sensors)
    else:
        if num_function_sensors > 0:
            h_function = 1. / num_function_sensors
            x_function = np.linspace(70, 90, num_function_sensors, endpoint=False) + 0.5 * h_function
        else:
            x_function = np.empty(shape=(0, 0))

        if num_derivative_sensors > 0:
            h_derivative = 1. / num_derivative_sensors
            x_derivative = np.linspace(0, 1, num_derivative_sensors, endpoint=False) + 0.5 * h_derivative
        else:
            x_derivative = np.empty(shape=(0, 0))

    all_a = np.array([280.+a for a in range(num_a)])
    data = {
        "all_a": all_a.tolist(),
        "b": b,
        "c": c,
        "x_function": x_function.tolist(),
        "sigma_noise_function": sigma_noise_function,
        "x_derivative": x_derivative.tolist(),
        "sigma_noise_derivative": sigma_noise_derivative,
        "seed": seed
    }
    with open(virtual_experimental_metadata_file, "w") as f:
        yaml.dump(data, f, default_flow_style=None)


def main():
    # create metadata for an exactly linear model (no model bias)
    write_metadata_to_yaml("virtual_experiment_linear_model_meta.yaml",
                           num_a=1, b=5, c=0,
                           num_function_sensors=10, sigma_noise_function=0.,
                           num_derivative_sensors=0, sigma_noise_derivative=0.,
                           center=False)

    write_metadata_to_yaml("virtual_experiment_linear_model_with_noise_meta.yaml",
                           num_a=1, b=5, c=0,
                           num_function_sensors=10, sigma_noise_function=0.01,    # 1.0
                           num_derivative_sensors=0, sigma_noise_derivative=0.01,
                           center=False)
'''
    # create metadata for a quadratic model (thus the linear model has a model bias)
    write_metadata_to_yaml("virtual_experiment_quadratic_model_meta.yaml",
                           num_a=1, b=5, c=1,
                           num_function_sensors=10, sigma_noise_function=0.,
                           num_derivative_sensors=0, sigma_noise_derivative=0.,
                           center=True)
'''

if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    main()
