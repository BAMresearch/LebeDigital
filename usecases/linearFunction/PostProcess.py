from Correlation import *
import pickle

def postprocess(trace, virtual_experiment, configurations, forward_model, data_base, out=None):
    s = pm.summary(trace)
    means = s["mean"]
    colors = cm.rainbow(np.linspace(0, 1, len(configurations)))
    for color, configuration in zip(colors,configurations):
        #join the two dictionaries
        all_prms = {**means.to_dict(),**configuration.parameters}
        model_response = forward_model(all_prms, configuration.sensor_groups)

        '''
        compute the posterior predictive using a sampling approach
        that is the response for a new x
        '''
        b_tr = trace.get_values('b', combine=True)
        noise_std_function_tr = trace.get_values('sigma_std_function', combine=True)
        noise_std_derivative_tr = trace.get_values('sigma_std_derivative', combine=True)
        # this is suboptimal, since the result has already been computed for all x (but is not stored in the trace)
        # since for the custom distribution pymc3 just needs the loglike, not the individual observations
        # could be fixed by internally serializing already computed model responses in a database and only
        # recompute if required
        # only if x is different than the training, this recomputation is necessary
        result_all_trace_with_noise = {}
        result_all_trace_without_noise = {}
        for b, noise_std_function, noise_std_derivative in zip(b_tr, noise_std_function_tr, noise_std_derivative_tr):
            single_result_without_noise = forward_model(
                {'a': configuration.parameters['a'], 'b': b}, configuration.sensor_groups)
            single_result_with_noise = {}
            for sensor_group in configuration.sensor_groups:
                #this is not fully correct since the noise is correlated, but for plotting, this should be fine
                if sensor_group.get_type()=="function":
                    single_result_with_noise[sensor_group] =  single_result_without_noise[sensor_group] + \
                        noise_std_function * np.random.randn(len(sensor_group.locations))
                if sensor_group.get_type()=="derivative":
                    single_result_with_noise[sensor_group] =  single_result_without_noise[sensor_group] + \
                        noise_std_derivative * np.random.randn(len(sensor_group.locations))

                if sensor_group in result_all_trace_without_noise:
                    result_all_trace_without_noise[sensor_group] = np.vstack([result_all_trace_without_noise[sensor_group],
                                                                              single_result_without_noise[sensor_group]])
                    result_all_trace_with_noise[sensor_group] = np.vstack(
                        [result_all_trace_with_noise[sensor_group],
                         single_result_with_noise[sensor_group]])
                else:
                    result_all_trace_without_noise[sensor_group] = single_result_without_noise[sensor_group]
                    result_all_trace_with_noise[sensor_group] = single_result_with_noise[sensor_group]

        for sensor_group in configuration.sensor_groups:
            if sensor_group.get_type()=="function" :
                plot = plt.figure(1,figsize=(5,2.5))
            if sensor_group.get_type()=="derivative" :
                plot = plt.figure(2,figsize=(5,2.5))
            #plt.title("sensor type = "+ sensor_group.get_type())
            if configuration.parameters['a'] == 4:
                plot = plt.fill_between(sensor_group.locations,
                                        np.percentile(result_all_trace_with_noise[sensor_group], 5, axis=0), \
                                        np.percentile(result_all_trace_with_noise[sensor_group], 95, axis=0), \
                                        color=color, alpha=0.06, label='90% Intervall')
            else:
                plot = plt.fill_between(sensor_group.locations,
                                        np.percentile(result_all_trace_with_noise[sensor_group], 5, axis=0), \
                                        np.percentile(result_all_trace_with_noise[sensor_group], 95, axis=0), \
                                        color=color, alpha=0.06)

            if configuration.parameters['a'] == 4:
                plot = plt.fill_between(sensor_group.locations,
                                        np.percentile(result_all_trace_without_noise[sensor_group], 5, axis=0), \
                                        np.percentile(result_all_trace_without_noise[sensor_group], 95, axis=0), \
                                        color=color, alpha=0.2, label='90% Intervall mit Noise')
            else:
                plot = plt.fill_between(sensor_group.locations,
                                        np.percentile(result_all_trace_without_noise[sensor_group], 5, axis=0), \
                                        np.percentile(result_all_trace_without_noise[sensor_group], 95, axis=0), \
                                        color=color, alpha=0.2)

            #estimated forward model with mean posterior parameters
            if configuration.parameters['a'] == 4:
                plt.plot(sensor_group.locations, model_response[sensor_group], linestyle='dashed', linewidth=1,
                         color=color,label = 'Mittelwert Prognose')
            else:
                plt.plot(sensor_group.locations, model_response[sensor_group], linestyle='dashed', linewidth=1,
                         color=color)
            #exact virtual experiment
            if configuration.parameters['a'] == 4:
                plt.plot(sensor_group.locations, virtual_experiment.result(sensor_group.get_type(), configuration,
                                                                           sensor_group.locations), color=color,
                                                                            label = 'Virtuelles Experiment')
            else:
                plt.plot(sensor_group.locations, virtual_experiment.result(sensor_group.get_type(), configuration,
                                                                           sensor_group.locations), color=color)

            #noisy virtual experiment
            if configuration.parameters['a'] == 4:
                plt.scatter(sensor_group.locations, data_base.data[configuration][sensor_group], color=color, label = 'Daten',s=1)
            else:
                plt.scatter(sensor_group.locations, data_base.data[configuration][sensor_group], color=color, s=1)
    plt.legend(fontsize='xx-small',bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig(out, dpi=600)
    plt.clf()
    plt.legend(fontsize='xx-small',bbox_to_anchor=(.0, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.hist(b_tr, bins=30, density=True)
    plt.savefig(os.path.splitext(out)[0]+'_b_post.png',dpi=600)
    plt.clf()
    plt.hist(noise_std_function_tr, bins=30, density=True)
    plt.savefig(os.path.splitext(out)[0]+'_noise_post.png',dpi=600)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import sys

    try:
        options_file = sys.argv[1]
    except:
        options_file = "NoCorrelation_coarse.yaml"

    trace_file = Path(options_file).with_suffix(".pkl")
    plot_file = Path(options_file).with_suffix(".png")

    [model_error, virtual_experiment, configurations, forward_model, data_base] = setup(Opts(options_file))

    with open(trace_file, 'rb') as buff:
        trace_data = pickle.load(buff)

    trace = trace_data['trace']

    postprocess(trace, virtual_experiment, configurations, forward_model, data_base, out=plot_file)
