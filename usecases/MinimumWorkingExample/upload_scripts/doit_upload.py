import logging
import os
import sys
from collections import defaultdict
from math import isnan
from pathlib import Path
import pandas as pd

import yaml

from lebedigital.openbis.interbis import Interbis


def upload_to_openbis_doit(
        metadata_path: str,
        processed_data_path: str,
        raw_data_path: str,
        mixture_metadata_file_path: str,
        mixture_data_path: str,
        output_path: str,
        config: dict):
    """Function for uploading data to the openbis datastore from within the doit environment

    Needed parameters in the config dict are:

    'datastore_url': Url to the openBIS datastore
    'space': Space within openbis for the sample
    'project': Project under specified space for the sample
    'emodul_collection': Collection under specified project for the emodul sample
    'mixture_collection': Collection under specified project for the mixture sample
    'sample_code': Code for the new type of the sample
    'sample_prefix': Prefix for the new type of the sample
    'mixture_code': Code for the new type of the mixture
    'mixture_prefix': Prefix for the new type of the mixture
    'verbose': If true the output will be printed to console, optional
    'runson': Specifies if the function is running on GitHub Actions or locally.
        To be parsed from command line by running 'doit runson=notactions'.

    Args:
        metadata_path (str): Path to the metadata yaml file
        processed_data_path (str): Path to the processed data file
        raw_data_path (str): Path to the raw data file
        mixture_metadata_file_path (str): Path to the mixture metadata file
        mixture_data_path (str): Path to the mixture data file
        output_path (str): Path where the samples overview should be saved
        config (dict): A dictionary containing the necessary info for uploading to openbis
    """

    # This does not want to work with the doit environment, some stuff gets printed while some does not
    if 'verbose' in config and config['verbose']:
        sys.stdout = sys.__stdout__
        logging.basicConfig(level=logging.DEBUG)
    else:
        sys.stdout = open(os.devnull, 'w')

    # Skipping the upload if the platform the code runs on can't log in to the Datastore
    if config['runson'] == 'actions':
        args = locals()
        _actions_run(
            metadata_path=args["metadata_path"],
            processed_data_path=args["processed_data_path"],
            raw_data_path=args["raw_data_path"],
            mixture_metadata_file_path=args["mixture_metadata_file_path"],
            mixture_data_path=args["mixture_data_path"],
            output_path=args["output_path"],
            config=args["config"]
        )
        return

    # Connecting to the datastore
    o = Interbis(config['datastore_url'])
    o.connect_to_datastore()

    # Setting "constants"
    _SPACE = config['space']
    _PROJECT = config['project']
    _EMODUL_COLLECTION = f"/{_SPACE}/{_PROJECT}/{config['emodul_collection']}"
    _MIXTURE_COLLECTION = f"/{_SPACE}/{_PROJECT}/{config['mixture_collection']}"

    """
    MIXTURE EXPERIMENTAL STEP FROM HERE ON
    """

    # We skip the mixture upload when the mixture yaml is not found
    if mixture_metadata_file_path:
        # Reading the metadata from output metadata yaml file
        mixture_metadata = _read_metadata_emodul(mixture_metadata_file_path)

        # Converting NaN values to 0.0 as openBIS does not accept NaNs
        for key, val in mixture_metadata.items():
            if type(val) == float and isnan(val):
                mixture_metadata[key] = 0.0

        # Transforming the metadata yaml file into format accepted by o.create_sample_type()
        mixture_metadata_types_dict = _reformat_sample_dict(mixture_metadata)

        mixture_sample_type = _create_mixture_sample_type(o, config=config,
                                                          sample_type_dict=mixture_metadata_types_dict)

        mixture_sample = _mixture_upload(
            o,
            mixture_metadata_dict=mixture_metadata,
            sample_name=os.path.splitext(os.path.basename(mixture_metadata_file_path))[0],
            mixture_sample_type=mixture_sample_type.code,
            mixture_data_filepath=mixture_data_path,
            space=_SPACE,
            project=_PROJECT,
            collection=_MIXTURE_COLLECTION
        )
    else:
        mixture_sample = "Not Found"

    """
    EMODUL EXPERIMENTAL STEP FROM HERE ON
    """

    # Reading the metadata from output metadata yaml file
    emodul_metadata = _read_metadata_emodul(metadata_path)
    for key, val in emodul_metadata.items():
        if type(val) == float and isnan(val):
            emodul_metadata[key] = 0.0

    # Converting NaN values to 0.0 as openBIS does not accept NaNs
    emodul_metadata_type_dict = _reformat_sample_dict(emodul_metadata)

    # Creating the emodul sample type with the formatted dict
    emodul_sample_type = o.create_sample_type(
        sample_code=config['sample_code'],
        sample_prefix=config['sample_prefix'],
        sample_properties=emodul_metadata_type_dict,
    )
    logging.debug(f'Emodul Sample Type {emodul_sample_type.code} created')

    # Initializing the new emodul sample
    emodul_sample = o.new_sample(
        type=emodul_sample_type.code,
        space=_SPACE,
        project=_PROJECT,
        collection=_EMODUL_COLLECTION,
        parents=mixture_sample.identifier
    )

    # Setting the metadata from the yaml file metadata, setting '$name' to 'experimentName'
    emodul_sample.set_props(emodul_metadata)
    emodul_sample.set_props({'$name': emodul_sample.props.all()['experimentname']})
    logging.debug(emodul_sample.metadata)

    # Checking if a sample with that name was already uploaded to that space in the datastore
    # Reason: Not uploading duplicate samples. Will check for duplicates based on '$name' property
    exist_emodul_sample_df = o.get_samples(
        space=_SPACE,
        project=_PROJECT,
        collection=_EMODUL_COLLECTION,
        props='$name').df

    exist_emodul_sample_df.columns = exist_emodul_sample_df.columns.str.lower()
    exist_emodul_sample_list = exist_emodul_sample_df['$name'].to_list()

    # Getting the name from props
    emodul_sample_name = emodul_sample.props.all()['$name']

    # If sample name found in all samples under that directory -> skip upload and fetch the sample instead
    # Else -> Upload the created sample to the Datastore
    if emodul_sample_name not in exist_emodul_sample_list:
        emodul_sample.save()
        logging.debug(f'Emodul Sample {emodul_sample.code} uploaded')
    else:
        # We checked that it exists, so we throw away the old sample and replace it with a fetched
        # sample from the datastore. The Identifier of the sample will be in column one row one of the dataframe
        emodul_sample = o.get_sample(exist_emodul_sample_df.iloc[0, 0])
        logging.debug(f'Emodul Sample {emodul_sample.code} found in datastore')

    # Naming the datasets
    raw_dataset_name = emodul_sample.props.all()['$name'] + '_raw'
    processed_dataset_name = emodul_sample.props.all()['$name'] + '_processed'

    # Checking whether the datasets already exist under the emodul sample
    emo_exist_dataset_df = o.get_datasets(sample=emodul_sample.identifier, props='$name').df
    emo_exist_dataset_df.columns = emo_exist_dataset_df.columns.str.lower()

    # Most sketchy stuff in pybis, I have no idea when the $NAME has to be capitalized and when it does not.
    # I suspect it depends on the object for which the get is called (sample/dataset/collection/... object)
    emo_exist_dataset_list = emo_exist_dataset_df['$name'].to_list()

    # Creating a transaction to upload the datasets at the same time

    # Raw Dataset Upload
    if raw_dataset_name not in emo_exist_dataset_list:
        raw_emo_dataset_props = {'$name': raw_dataset_name}

        emo_raw_dataset = o.new_dataset(
            type='RAW_DATA',
            collection=_EMODUL_COLLECTION,
            sample=emodul_sample.identifier,
            files=[raw_data_path],
            props=raw_emo_dataset_props
        )
        # Adding to transaction
        emo_raw_dataset.save()
        logging.debug('Raw Emodul Dataset uploaded')
    else:
        # The _ can be changed to a variable if the dataset is needed
        # _ = o.get_dataset(permIds=emo_exist_dataset_df.iloc[0, 0])
        logging.debug('Raw Emodul Dataset fetched from datastore')

    # Processed Dataset Upload
    if processed_dataset_name not in emo_exist_dataset_list:
        processed_emo_dataset_props = {'$name': processed_dataset_name}

        emo_processed_dataset = o.new_dataset(
            type='PROCESSED_DATA',
            collection=_EMODUL_COLLECTION,
            sample=emodul_sample.identifier,
            files=[processed_data_path],
            props=processed_emo_dataset_props
        )
        # Adding to transaction
        emo_processed_dataset.save()
        logging.debug('Processed Emodul Dataset uploaded')
    else:
        # The _ can be changed to a variable if the dataset is needed
        # _ = o.get_dataset(permIds=emo_exist_dataset_df.iloc[0, 0])
        logging.debug('Processed Emodul Dataset fetched from datastore')

    # We load the object from the datastore before printing as a sort of manual check
    # if the function worked as it was supposed to.
    emo_output_sample = o.get_sample(emodul_sample.identifier)

    file_name_with_extension = emo_output_sample.code + '.yaml'
    logfile_path = str(Path(output_path, file_name_with_extension))

    # Saving the log of the sample to a file in the output for dodo to have something to check
    with open(logfile_path, 'w') as file:
        print(emo_output_sample, file=file)

    sys.stdout = sys.__stdout__


def _read_metadata_emodul(yaml_path: str):
    """Reads the metadata as it is saved in the emodul_metadata directory

    Args:
        yaml_path (str): path to the yaml file
    """

    with open(yaml_path, 'r') as file:
        data = yaml.safe_load(file)
        data = {k.lower(): v for k, v in data.items()}
        return data


def _reformat_sample_dict(loaded_dict: dict):
    conv_dict = {
        str: 'VARCHAR',
        float: 'REAL',
        int: 'INTEGER',
    }

    output_dict = defaultdict(lambda: ["NA", "NA", "NA"])
    output_dict['$name'] = ['VARCHAR', 'Name', 'Name']

    for key, val in loaded_dict.items():
        output_dict[key.lower()] = [conv_dict[type(val)], f"{key}_label", f"{key}_description"]

    return dict(output_dict)


def _create_mixture_sample_type(o: Interbis, config: dict, sample_type_dict: dict):
    # Creating the mixture sample type with the formatted dict
    mixture_sample_type = o.create_sample_type(
        sample_code=config['mixture_code'],
        sample_prefix=config['mixture_prefix'],
        sample_properties=sample_type_dict,
    )
    logging.debug(f'mixture type created: {mixture_sample_type.code}')
    return mixture_sample_type


def _mixture_upload(
        o: Interbis,
        mixture_metadata_dict: dict,
        sample_name: str,
        mixture_sample_type: str,
        mixture_data_filepath: str,
        space: str,
        project: str,
        collection: str):
    # Initializing the new mixture sample
    mixture_sample = o.new_sample(
        type=mixture_sample_type,
        space=space,
        project=project,
        collection=collection
    )

    # Setting the props from metadata and adding '$name' for better readability in the web view
    mixture_sample.set_props(mixture_metadata_dict)

    # mixture_sample.set_props({'$name': os.path.splitext(os.path.basename(mixture_metadata_file_path))[0]})
    mixture_sample.set_props({"$name": sample_name})
    logging.debug(mixture_sample.props)

    # Checking if a sample with that name was already uploaded to that space in the datastore
    # Reason: Not uploading duplicate samples. Will check for duplicates based on '$name' property
    exist_mixture_sample_df = o.get_samples(
        space=space,
        project=project,
        collection=collection,
        props='$name').df

    exist_mixture_sample_df.columns = exist_mixture_sample_df.columns.str.lower()
    exist_mixture_sample_list = exist_mixture_sample_df['$name'].to_list()

    # Getting the name from props
    mixture_sample_name = mixture_sample.props.all()['$name']

    # If sample name found in all samples under that directory -> skip upload and fetch the sample instead
    # Else -> Upload the created sample to the Datastore
    if mixture_sample_name not in exist_mixture_sample_list:
        mixture_sample.save()
        logging.debug(f'mixture {mixture_sample.code} uploaded')
    else:
        # We checked that it exists, so we throw away the old sample and replace it with a fetched
        # sample from the datastore. The Identifier of the sample will be in column one row one of the dataframe
        mixture_sample = o.get_sample(exist_mixture_sample_df.iloc[0, 0])
        logging.debug(f'mixture found in dataset')

    mixture_dataset_name = mixture_sample_name + '_raw'

    # Same as before, checking whether the dataset was uploaded already
    mix_exist_dataset_df = o.get_datasets(sample=mixture_sample.identifier, props=['$name']).df
    mix_exist_dataset_df.columns = mix_exist_dataset_df.columns.str.lower()

    # Most sketchy stuff in pybis, I have no idea when the $NAME has to be capitalized and when it does not.
    # I suspect it depends on the object for which the get is called (sample/dataset/collection/... object)
    mix_exist_dataset_list = mix_exist_dataset_df['$name'].to_list()

    if mixture_dataset_name not in mix_exist_dataset_list:
        raw_mix_dataset_props = {'$name': mixture_dataset_name}

        mixture_dataset = o.new_dataset(
            type='RAW_DATA',
            collection=collection,
            sample=mixture_sample.identifier,
            files=[mixture_data_filepath],
            props=raw_mix_dataset_props
        )
        mixture_dataset.save()
        logging.debug('Mixture Dataset uploaded')
    else:
        # Can change the _ to some variable if the dataset becomes needed
        # _ = o.get_dataset(permIds=mix_exist_dataset_df.iloc[0, 0])
        logging.debug('Mixture Dataset fetched from datastore')

    return mixture_sample


def _actions_run(
        metadata_path: str,
        processed_data_path: str,
        raw_data_path: str,
        mixture_metadata_file_path: str,
        mixture_data_path: str,
        output_path: str,
        config: dict):
    logging.debug("Running in no database connection mode")

    mix_file_yaml = os.path.basename(mixture_metadata_file_path) if mixture_metadata_file_path else "Not Found"
    mix_file_data = os.path.basename(mixture_data_path) if mixture_metadata_file_path else "Not Found"
    emodul_file_yaml = os.path.basename(metadata_path)
    emodul_raw_data = os.path.basename(raw_data_path)
    emodul_processed_data = os.path.basename(processed_data_path)

    logging.debug(mix_file_yaml)
    logging.debug(mix_file_data)
    logging.debug(emodul_file_yaml)
    logging.debug(emodul_raw_data)
    logging.debug(emodul_processed_data)

    output_dict = {
        'ran_on': 'github_actions',
        'db_url': config["datastore_url"],
        'mix_file_yaml': mix_file_yaml,
        'mix_file_data': mix_file_data,
        'emodul_file_yaml': emodul_file_yaml,
        'emodul_raw_data': emodul_raw_data,
        'emodul_processed_data': emodul_processed_data
    }
    file_name_with_extension = str(os.path.splitext(os.path.basename(emodul_file_yaml))[0] + '_oB_upload_log.yaml')
    with open(Path(output_path, file_name_with_extension), 'w') as file:
        _ = yaml.dump(output_dict, file)
    return