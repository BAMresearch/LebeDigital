from lebedigital.openbis.expstep import ExpStep
from lebedigital.openbis.interbis import Interbis
from pathlib import Path
import sys
import os
import yaml

def upload_to_openbis_doit(metadata_path: str, processed_data_path: str, raw_data_path: str, output_path: str, config: dict):
    """Function for uploading data to the openbis datastore from within te doit environment

    Needed parameters in the config dict are:

    'datastore_url': Url to the openBIS datastore
    'space': Space within openbis for the sample
    'project': Project under specified space for the sample
    'collection': Collection under specified project for the sample
    'sample_code': Code for the new type of the sample
    'sample_prefix': Prefix for the new type of the sample
    'verbose': If true the output will be printed to console, optional
    'runson': Specifies if the function is running on github actions or locally. To be parsed from command line by running 'doit runson=notactions'.

    Args:
        metadata_path (str): Path to the metadata yaml file
        processed_data_path (str): Path to the processed data file
        raw_data_path (str): Path to the raw data file
        output_path (str): Path where the samples overview should be saved
        config (dict): A dictionary containing the necessary info for uploading to openbis
    """

    # This does not want to work with the doit environment, some stuff gets printed while some does not
    if 'verbose' in config and config['verbose']:
        sys.stdout = sys.__stdout__
    else:
        sys.stdout = open(os.devnull, 'w')

    if config['runson'] == 'actions':
        output_dict = {'ran_on': 'github_actions'}
        file_name_with_extension = 'logfile.yaml'
        with open(Path(output_path, file_name_with_extension), 'w') as file:
            documents = yaml.dump(output_dict, file)
        return

    emodul_sample = ExpStep(
        name=os.path.splitext(os.path.basename(metadata_path))[0],
        space=config['space'],
        project=config['project'],
    )

    o = Interbis(config['datastore_url'])
    o.connect_to_datastore()

    emodul_sample.collection = emodul_sample.find_collection(
        o,
        config['collection'],
        id_type=1,
    )

    emodul_sample.sync_name(get_from='name')

    emodul_sample.read_metadata_emodul(metadata_path)

    emodul_sample_type = ExpStep.create_sample_type_emodul(
        o,
        sample_code=config['sample_code'],
        sample_prefix=config['sample_prefix'],
        sample_properties=emodul_sample.metadata,
    )

    emodul_sample.type = emodul_sample_type.code

    emodul_sample.upload_expstep(o)

    processed_dataset_name = emodul_sample.name + '_processed'

    emodul_sample.upload_dataset(
        o,
        props={
            '$name': processed_dataset_name,
            'files': [processed_data_path],
            'data_type': 'PROCESSED_DATA'
        }
    )

    raw_dataset_name = emodul_sample.name + '_raw'

    emodul_sample.upload_dataset(
        o,
        props={
            '$name': raw_dataset_name,
            'files': [raw_data_path],
            'data_type': 'RAW_DATA'
        }
    )

    # We load the object from the datastore before printing as a sort of manual check if the function worked as it was supposed to.
    output_sample = ExpStep.load_sample(o, emodul_sample.identifier)
    # output_sample.info()

    file_name_with_extension = output_sample.name + '.yaml'

    output_sample.save_sample_yaml(Path(output_path, file_name_with_extension))

    sys.stdout = sys.__stdout__