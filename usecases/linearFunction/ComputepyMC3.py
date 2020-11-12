from Correlation import *
import pickle

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import sys

    try:
        options_file = sys.argv[1]
    except:
        options_file = "NoCorrelation_coarse.yaml"

    pickle_file = Path(options_file).with_suffix(".pkl")

    [model_error, virtual_experiment, configurations, forward_model, data_base] = setup(Opts(options_file))
    trace = run_pymc3(model_error)

    with open(pickle_file, 'wb') as buff:
        pickle.dump({'trace': trace}, buff)
