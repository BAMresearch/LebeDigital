import os
from pathlib import Path

import yaml
from doit import create_after
from doit.task import clean_targets

from lebedigital.raw_data_processing.youngs_modulus_data \
    .emodul_metadata_extraction import emodul_metadata

from lebedigital.mapping.emodul_mapping import generate_knowledge_graph

from lebedigital.raw_data_processing.youngs_modulus_data \
    .emodul_generate_processed_data import processed_data_from_rawdata

from lebedigital.raw_data_processing.mixture \
    .mixture_metadata_extraction import extract_metadata_mixture

from upload_scripts.doit_upload import upload_to_openbis_doit

from doit import get_var

# set a variable to define a cheap or full run
# the default "doit" is set to "doit mode=cheap"
# "cheap" option is to reduce computation time once the workflow gets more expensive (calibration)
#  - currently this means: all mixes, one tests data + KG
# "single" option is to test the dodo file on a single example (similar to cheap but only a single mix)
# any other mode value runs the expensive version i.e. "doit mode=full"
config = {"mode": get_var('mode', 'cheap')}

# when "cheap option" or "single" is run, only this souce of raw data is processed
single_example_name = 'Wolf 8.2 Probe 1'
# TODO: (if we want to keep using a single example) automatic identification of corresponding mix
single_mix_name = '2014_12_10 Wolf.xls'  # corresponding mix for the "single" example

DOIT_CONFIG = {'verbosity': 2}

# openbis config needed for the upload to the datastore
openbis_config = {
    'datastore_url': get_var("url", 'https://test.datastore.bam.de/openbis/'),
    'space': 'CKUJATH',
    'project': 'LEBEDIGITAL',
    'emodul_collection': 'LEBEDIGITAL_EMODUL_COLLECTION',
    'mixture_collection': 'LEBEDIGITAL_MIXTURE_COLLECTION',
    'sample_code': 'EXPERIMENTAL_STEP_EMODUL',
    'sample_prefix': 'EMODUL',
    'mixture_code': 'EXPERIMENTAL_STEP_EMODUL_MIX',
    'mixture_prefix': 'EMODUL_MIX',
    'verbose': False,
    # if actions is specified the task will be completed but the openbis connection will be skipped
    # we need to skip the openbis functions on github actions as they need a password to run
    'runson': get_var('runson', 'actions'),
}

# parent directory of the minimum working example
ParentDir = os.path.dirname(Path(__file__))

# EMODULE PATHS 
# defining paths for emodule
emodul_output_directory = Path(ParentDir, 'emodul')  # folder with metadata yaml files
raw_data_emodulus_directory = Path(ParentDir, 'Data', 'E-modul')  # folder with folders of raw data files
metadata_emodulus_directory = Path(emodul_output_directory, 'metadata_yaml_files')  # folder with metadata yaml files
processed_data_emodulus_directory = Path(emodul_output_directory, 'processed_data')  # folder with csv data files
knowledge_graphs_directory = Path(emodul_output_directory, 'knowledge_graphs')  # folder with KG ttl files
openbis_samples_directory = Path(emodul_output_directory, 'openbis_samples')

# create folder, if it is not there
Path(emodul_output_directory).mkdir(parents=True, exist_ok=True)

# MIXTURE PATHS
# defining paths for mixture
raw_data_mixture_directory = Path(ParentDir, 'Data', 'Mischungen')  # folder with raw data files (excel)
mixture_output_directory = Path(ParentDir, 'mixture')  # folder with folders
metadata_mixture_directory = Path(mixture_output_directory,
                                  'metadata_yaml_files')  # folder with mixture metadata yaml files
mixture_knowledge_graphs_directory = Path(mixture_output_directory, 'knowledge_graphs')  # folder with KG ttl files

# List with mixes with problems, to be excluded for now
excluded_mix_list = ['2014_08_04 Rezepturen_auf 85 Liter_Werner_Losert.xlsx']

# create folder, if it is not there
Path(mixture_output_directory).mkdir(parents=True, exist_ok=True)


# TASKS
# extract metadata for the mixture
def task_extract_metadata_mixture():
    # create folder, if it is not there
    Path(metadata_mixture_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'single':
        list_raw_data_mixture_files = [Path(raw_data_mixture_directory, single_mix_name)]

    else:  # make a list of all files
        list_raw_data_mixture_files = os.scandir(raw_data_mixture_directory)

    for f in list_raw_data_mixture_files:
        if f.is_file():
            raw_data_path = Path(f)
            yaml_metadata_file = Path(metadata_mixture_directory, f.name.split(".xls")[0] + '.yaml')

            if f.name not in excluded_mix_list:
                yield {
                    'name': f.name,
                    'actions': [(extract_metadata_mixture, [raw_data_path, metadata_mixture_directory])],
                    'targets': [yaml_metadata_file],
                    'clean': [clean_targets]
                }


# extract standardized meta data for Young' modulus tests
def task_extract_metadata_emodul():
    # create folder, if it is not there
    Path(metadata_emodulus_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap' or config['mode'] == 'single':
        list_raw_data_emodulus_directories = [Path(raw_data_emodulus_directory, single_example_name)]
    else:  # go through all files
        list_raw_data_emodulus_directories = os.scandir(raw_data_emodulus_directory)

    for f in list_raw_data_emodulus_directories:
        if f.is_dir():
            raw_data_path = Path(f)
            raw_data_file = Path(f, 'specimen.dat')
            yaml_metadata_file = Path(
                metadata_emodulus_directory, f.name + '.yaml')
            yield {
                'name': f.name,
                'actions': [(emodul_metadata, [raw_data_path, yaml_metadata_file])],
                'file_dep': [raw_data_file],
                'targets': [yaml_metadata_file],
                'clean': [clean_targets]
            }


# extract standardized processed data for Young' modulus tests


def task_extract_processed_data_emodul():
    # create folder, if it is not there
    Path(processed_data_emodulus_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap' or config['mode'] == 'single':
        list_raw_data_emodulus_directories = [Path(raw_data_emodulus_directory, single_example_name)]
    else:  # go through all files
        list_raw_data_emodulus_directories = os.scandir(raw_data_emodulus_directory)

    for f in list_raw_data_emodulus_directories:
        if f.is_dir():
            raw_data_file = Path(f, 'specimen.dat')
            # the name of the csv file is the file name of the raw data
            # is processed_data_directory + directory_raw_data.csv
            csv_data_file = Path(
                processed_data_emodulus_directory, f.name + '.csv')

            yield {
                'name': f.name,
                'actions': [(processed_data_from_rawdata, [f, csv_data_file])],
                'file_dep': [raw_data_file],
                'targets': [csv_data_file],
                'clean': [clean_targets]
            }


# generate knowledgeGraphs
@create_after(executed='extract_metadata_emodul')
def task_export_knowledgeGraph_emodul():
    # create folder, if it is not there
    Path(knowledge_graphs_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap' or config['mode'] == 'single':
        list_metadata_yaml_files = [Path(metadata_emodulus_directory, single_example_name + '.yaml')]
    else:  # go through all files
        # list of all meta data files....
        list_metadata_yaml_files = os.scandir(metadata_emodulus_directory)

    # check directory, if

    for f in list_metadata_yaml_files:
        if f.is_file():
            # path to metadata yaml
            metadata_file_path = Path(f)
            name_of_ttl = f.name.replace('.yaml', '.ttl')
            name_of_cvs = f.name.replace('.yaml', '.csv')
            # path the processed data csv
            processed_data_file_path = Path(
                processed_data_emodulus_directory, name_of_cvs)
            # path to output file KG
            knowledge_graph_file = Path(
                knowledge_graphs_directory, name_of_ttl)

            yield {
                'name': name_of_cvs,
                'actions': [(generate_knowledge_graph, [metadata_file_path,
                                                        knowledge_graph_file])],
                'file_dep': [metadata_file_path, processed_data_file_path],
                'targets': [knowledge_graph_file],
                'clean': [clean_targets]
            }


@create_after(target_regex='.*emodul$')
def task_upload_to_openbis():
    # create folder, if it is not there
    Path(openbis_samples_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap' or config['mode'] == 'single':
        list_metadata_emodul_yaml_files = [Path(metadata_emodulus_directory, single_example_name)]
    else:  # go through all files
        list_metadata_emodul_yaml_files = os.scandir(metadata_emodulus_directory)

    for meta_f in list_metadata_emodul_yaml_files:
        # getting path for metadata yaml file
        metadata_file_path = Path(meta_f)

        # getting the processed file path
        processed_file_name = str(os.path.splitext(os.path.basename(metadata_file_path))[0])+'.csv'
        processed_file_path = Path(processed_data_emodulus_directory, processed_file_name)

        # here we get the corresponding mix_file
        with open(metadata_file_path, 'r') as file:
            data = yaml.safe_load(file)
            mix_file = data['mix_file']

        print(mix_file)
        # If mix file is found, find it, else set the paths to empty strings values indicating missing data
        if mix_file:
            mixture_metadata_file_path = Path(mixture_output_directory,
                                              'metadata_yaml_files',
                                              str(os.path.splitext(os.path.basename(mix_file))[0])+'.yaml')

            # raw_data_mixture_directory
            mixture_data_path = Path(raw_data_mixture_directory,
                                     str(os.path.splitext(os.path.basename(mix_file))[0]) + '.xls')
        else:
            mixture_metadata_file_path = ""
            mixture_data_path = ""

        # the raw data file is the specimen file in the corresponding folder in raw_data_emodulus_directory
        raw_data_file = Path(
            raw_data_emodulus_directory,
            os.path.splitext(os.path.basename(meta_f))[0], 'specimen.dat')

        sample_file_name = str(os.path.splitext(os.path.basename(metadata_file_path))[0]+'_oB_upload_log.yaml')
        sample_file_path = Path(openbis_samples_directory, sample_file_name)

        yield {
            'name': metadata_file_path,
            'actions': [(upload_to_openbis_doit,
                         [metadata_file_path, processed_file_path, raw_data_file,
                          mixture_metadata_file_path, mixture_data_path, openbis_samples_directory, openbis_config])],
            'file_dep': [metadata_file_path, processed_file_path],
            'targets': [sample_file_path],
            'clean': [clean_targets],
        }
