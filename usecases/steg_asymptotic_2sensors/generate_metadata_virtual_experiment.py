import numpy as np
import yaml


def write_metadata_to_yaml(virtual_experimental_metadata_file,
                           a, b, c, d, e, mA, mB,
                           num_Mat_sensors, sigma_noise_Mat,
                           center=False, seed=42):
    """Create meta data of virtual experiment
    Args:
        virtual_experimental_metadata_file: filename for metadata to be generated
        a: MatB coefficient of virtual experiment
        b: MatB coefficient of virtual experiment
        c: MatC coefficient of virtual experiment
        d: MatC coefficient of virtual experiment
        e: MatC coefficient of virtual experiment
        mA: model coefficient of Mat_B_C asymptotic model (=MatB/MatC) mA + mB/x
        mB: model coefficient of Mat_B_C asymptotic model (=MatB/MatC) mA + mB/x
        num_Mat_sensors(int): number of Mat sensors in the interval [x_start,x_end]
        sigma_noise_Mat: std deviation of the MatB noise

        center: if False, include start and end point, if True, divide in equal intervals and use midpoint
        seed: seed for the random number generator
    """

    x_start = 300
    x_end = 750

    if center is False:
        x_Mat_sensors = np.linspace(x_start, x_end, num_Mat_sensors)
    else:
        if num_Mat_sensors > 0:
            h_Mat = 1. / num_Mat_sensors
            x_Mat_sensors = np.linspace(x_start, x_end, num_Mat_sensors, endpoint=False) + 0.5 * h_Mat
        else:
            x_Mat_sensors = np.empty(shape=(0, 0))


#    all_a = np.array([1+2*a for a in range(num_a)])
    data = {
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "e": e,
        "mA": mA,
        "mB": mB,
        "x_Mat_sensors": x_Mat_sensors.tolist(),
        "sigma_noise_Mat": sigma_noise_Mat,
        "seed": seed
    }
    with open(virtual_experimental_metadata_file, "w") as f:
        yaml.dump(data, f, default_flow_style=None)


def main():
    # create metadata for an exactly asymptotic model (no model bias)
    write_metadata_to_yaml("virtual_experiment_asymptotic_model_meta.yaml",
                           a=-280, b=0.99, c=.0007, d=1., e=0.1, mA=0.99, mB=-280,
                           num_Mat_sensors=10, sigma_noise_Mat=0.,
                           center=False) # 20 num sensors

    write_metadata_to_yaml("virtual_experiment_asymptotic_model_with_noise_meta.yaml",
                           a=-280, b=0.99, c=.0007, d=1., e=0.1, mA=0.99, mB=-280,
                           num_Mat_sensors=10, sigma_noise_Mat=0.01,
                           center=False)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    main()
