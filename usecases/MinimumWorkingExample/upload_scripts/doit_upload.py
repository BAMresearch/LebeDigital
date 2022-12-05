from lebedigital.openbis.expstep import ExpStep
from lebedigital.openbis.interbis import Interbis
from pathlib import Path
import sys
import os
import yaml
import logging
from math import isnan


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

    if config['runson'] == 'actions':
        output_dict = {'ran_on': 'github_actions'}
        file_name_with_extension = 'logfile.yaml'
        with open(Path(output_path, file_name_with_extension), 'w') as file:
            documents = yaml.dump(output_dict, file)
        return

    o = Interbis(config['datastore_url'])
    o.connect_to_datastore()

    """
    MIXTURE EXPERIMENTAL STEP FROM HERE ON
    """

    mixture_sample = ExpStep(
        name=os.path.splitext(os.path.basename(mixture_metadata_file_path))[0],
        space=config['space'],
        project=config['project'],
    )

    mixture_sample.collection = o.get_collection_identifier(
        config['mixture_collection'],
    )

    print(mixture_metadata_file_path)
    mixture_sample.metadata = read_metadata_emodul(mixture_metadata_file_path)

    mixture_sample.metadata['$name'] = os.path.splitext(os.path.basename(mixture_metadata_file_path))[0]
    logging.debug(mixture_sample.metadata)

    # Creating the sample dict with data types (float, str, etc.) instead of the metadata values to create the
    # Sample Type for the mixture

    mixture_metadata_types_dict = {key: type(val) for key, val in mixture_sample.metadata.items()}

    mixture_sample_type = o.create_sample_type(
        sample_code=config['mixture_code'],
        sample_prefix=config['mixture_prefix'],
        sample_properties=mixture_metadata_types_dict,
    )

    mixture_sample.type = mixture_sample_type.code
    logging.debug(f'mixture type created: {mixture_sample.type}')

    for key, val in mixture_sample.metadata.items():
        if type(val) == float and isnan(val):
            mixture_sample.metadata[key] = 0.0

    mixture_sample.upload_expstep(o)
    logging.debug(f'mixture uploaded')

    mixture_dataset_name = mixture_sample.name + '_raw'

    mixture_sample.upload_dataset(
        o,
        props={
            '$name': mixture_dataset_name,
            'files': [mixture_data_path],
            'data_type': 'RAW_DATA'
        }
    )
    logging.debug('Mixture Dataset Uploaded')

    """
    EMODUL EXPERIMENTAL STEP FROM HERE ON
    """
    emodul_sample = ExpStep(
        name=os.path.splitext(os.path.basename(metadata_path))[0],
        space=config['space'],
        project=config['project'],
    )

    emodul_sample.collection = o.get_collection_identifier(
        config['emodul_collection'],
    )

    emodul_sample.metadata = read_metadata_emodul(metadata_path)

    emodul_sample.metadata['$name'] = emodul_sample.metadata['experimentname']
    logging.debug(emodul_sample.metadata)

    # Creating the sample dict with data types (float, str, etc.) instead of the metadata values to create the
    # Sample Type for the emodul sample

    emodul_sample_mixture_dict = {key: type(val) for key, val in emodul_sample.metadata.items()}

    emodul_sample_type = o.create_sample_type(
        sample_code=config['sample_code'],
        sample_prefix=config['sample_prefix'],
        sample_properties=emodul_sample_mixture_dict,
    )

    emodul_sample.type = emodul_sample_type.code
    logging.debug(f'sample type created: {emodul_sample.type}')

    # SETTING THE MIXTURE EXPSTEP AS A PARENT OF EMODUL EXPSTEP
    emodul_sample.parents = [o.get_sample_identifier(mixture_sample.name)]

    emodul_sample.upload_expstep(o)
    logging.debug('Sample Uploaded')

    processed_dataset_name = emodul_sample.name + '_processed'

    emodul_sample.upload_dataset(
        o,
        props={
            '$name': processed_dataset_name,
            'files': [processed_data_path],
            'data_type': 'PROCESSED_DATA'
        }
    )
    logging.debug('Processed Dataset Uploaded')

    raw_dataset_name = emodul_sample.name + '_raw'

    emodul_sample.upload_dataset(
        o,
        props={
            '$name': raw_dataset_name,
            'files': [raw_data_path],
            'data_type': 'RAW_DATA'
        }
    )
    logging.debug('Raw Dataset Uploaded')

    # We load the object from the datastore before printing as a sort of manual check
    # if the function worked as it was supposed to.

    output_sample = ExpStep.load_sample(o, emodul_sample.identifier)
    # output_sample.info()

    file_name_with_extension = output_sample.name + '.yaml'

    output_sample.save_sample_yaml(str(Path(output_path, file_name_with_extension)))

    sys.stdout = sys.__stdout__


def read_metadata_emodul(yaml_path: str):
    """Reads the metadata as it is saved in the emodul_metadata directory

    Args:
        yaml_path (str): path to the yaml file
    """

    with open(yaml_path, 'r') as file:
        data = yaml.safe_load(file)

        data = {k.lower(): v for k, v in data.items()}

        return data
