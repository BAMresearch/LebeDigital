import numpy as np
import yaml


def write_metadata_to_yaml(virtual_experimental_metadata_file,
                           num_a, b, c, num_function_sensors, num_derivative_sensors, center=False):
    """Create meta data of virtual experiment
    Args:
        virtual_experimental_metadata_file: filename for metadata to be generated
        num_a: number of constant offset a computed by 1+2*i with i in range(num_a)
        b: linear coefficient of the model
        c: quadratic coefficient of the model
        num_function_sensors(int): number of derivative sensors in the interval [0,1]
        num_derivative_sensors(int): number of derivative sensors in the interval [0,1]
        center: if False, include start and end point, if True, divide in equal intervals and use midpoint
    """

    if center is False:
        x_function = np.linspace(0, 1, num_function_sensors)
        x_derivative = np.linspace(0, 1, num_derivative_sensors)
    else:
        if num_function_sensors > 0:
            h_function = 1. / num_function_sensors
            x_function = np.linspace(0, 1, num_function_sensors, endpoint=False) + 0.5 * h_function
        else:
            x_function = np.empty(shape=(0, 0))

        if num_derivative_sensors > 0:
            h_derivative = 1. / num_derivative_sensors
            x_derivative = np.linspace(0, 1, num_derivative_sensors, endpoint=False) + 0.5 * h_derivative
        else:
            x_derivative = np.empty(shape=(0, 0))

    all_a = np.array([1+2*a for a in range(num_a)])
    data = {
        "all_a": all_a.tolist(),
        "b": b,
        "c": c,
        "x_function": x_function.tolist(),
        "x_derivative": x_derivative.tolist()
    }
    with open(virtual_experimental_metadata_file, "w") as f:
        yaml.dump(data, f, default_flow_style=None)


def main():
    # create metadata for an exactly linear model (no model bias)
    write_metadata_to_yaml("virtual_experiment_linear_model_meta.yaml",
        num_a=10, b=3, c=0, num_function_sensors=10, num_derivative_sensors=6, center=False)

    # create metadata for a quadratic model (thus the linear model has a model bias)
    write_metadata_to_yaml("virtual_experiment_quadratic_model_meta.yaml",
        num_a=1, b=0, c=2, num_function_sensors=1000, num_derivative_sensors=0, center=True)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    main()
