import os
import random
from pathlib import Path
import numpy as np
import pandas as pd

import yaml
import json
from doit import create_after, get_var
from doit.task import clean_targets
#from upload_scripts.doit_create_types import create_required_sample_types
#from upload_scripts.doit_upload import upload_to_openbis_doit

from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import \
    mix_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_generate_processed_data import \
    processed_data_from_rawdata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import \
    emodul_metadata
#from lebedigital.raw_data_processing.Compressive_strength.ComSt_generate_processed_data import \
    #processed_rawdata
#from lebedigital.raw_data_processing.Compressive_strength.ComSt_metadata_extraction import \
    #ComSt_metadata
from lebedigital.mapping.mappingscript import mapping

# from lebedigital.shacl import validation as shacl

# set a variable to define a cheap or full run
# the default "doit" is set to "doit mode=cheap"
# "cheap" option is to reduce computation time once the workflow gets more expensive (calibration)
#  - currently this means: all mixes, one tests data + KG
# "single" option is to test the dodo file on a single example (similar to cheap but only a single mix)
# any other mode value runs the expensive version i.e. "doit mode=full"
config = {"mode": get_var('mode', 'cheap')}

# when "cheap option" or "single" is run, only this source of raw data is processed
single_example_name = 'Wolf 8.2 Probe 1'
# TODO: (if we want to keep using a single example) automatic identification of corresponding mix
# corresponding mix sor the "single" example
single_mix_name = '20240220_7188_M01.xls'

single_KG_template = 'MixtureDesign_KG_Template.ttl'

DOIT_CONFIG = {'verbosity': 2}


# parent directory of the minimum working example
ParentDir = os.path.dirname(Path(__file__))


# EMODULE PATHS
# defining paths for emodule
# folder with metadat json files
emodul_output_directory = Path(ParentDir, 'emodul')
# folder with folders of raw data files
raw_data_emodulus_directory = Path(ParentDir, 'Data', 'E-modul')
# folder with metadata json files
metadata_emodulus_directory = Path(
    emodul_output_directory, 'mixture/metadata_json_files')
processed_data_emodulus_directory = Path(
    emodul_output_directory, 'processed_data')  # folder with csv data files
knowledge_graphs_output_directory = Path(
    ParentDir, 'Mapping_Example')  # folder with KG ttl files


# create folder, if it is not there
Path(emodul_output_directory).mkdir(parents=True, exist_ok=True)

# MIXTURE PATHS
# defining paths for mixture
# folder with raw data files (excel)
raw_data_mixture_directory = Path(ParentDir, 'Data', 'Mischungen_BAM')
mixture_output_directory = Path(ParentDir, 'mixture')   # folder with folders
metadata_mixture_directory = Path(mixture_output_directory, 'metadata_json_files')

# List with mixes with problems, to be excluded for now
excluded_mix_list = ['2014_08_04 Rezepturen_auf 85 Liter_Werner_Losert.xlsx']

# create folder, if it is not there
Path(mixture_output_directory).mkdir(parents=True, exist_ok=True)

# shacl directory where shapes are stored, to be tested against the created KGs
root_directory = Path(ParentDir).parent.parent
lebedigital_directory = Path(root_directory, 'lebedigital')
knowledge_graphs_Template_directory = Path(root_directory, 'lebedigital', 'ConcreteOntology')
shacl_directory = Path(lebedigital_directory, 'shacl')
# create folder, if it is not there
#union_output_path = Path(mixture_metadata_directory, "mixture_union_files")
#union_output_path.mkdir(parents=True, exist_ok=True)


def fixnames(path):
    for f in os.listdir(path):
        r = f.replace(" ","")
        if( r != f):
            os.rename(f,r)

# TASKS
# extract metadata for the mixture
def task_extract_metadata_mixture():
    # create folder, if it is not there
    Path(metadata_mixture_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'single':
        list_raw_data_mixture_files = [
            Path(raw_data_mixture_directory, single_mix_name)]

    else:  # make a list of all files
        list_raw_data_mixture_files = os.scandir(raw_data_mixture_directory)

    for f in list_raw_data_mixture_files:
        if f.is_file():
            raw_data_path = Path(f)
            # Extract the file name without extension and remove spaces
            file_name_without_extension = f.name.split(".xls")[0]
            #json_metadata_file = Path(metadata_mixture_directory, " ")
            json_metadata_file = str(metadata_mixture_directory) + "/"
            json_metadata = Path(json_metadata_file,  file_name_without_extension + '.json')
            if f.name not in excluded_mix_list:
                yield {
                    'name': f.name,
                    'actions': [(mix_metadata, [raw_data_path, json_metadata_file])],
                    'targets': [json_metadata],
                    'clean': [clean_targets]
                }
# extract standardized metadata for Young's modulus tests
#@create_after('extract_metadata_mixture')
#def task_extract_metadata_emodul():
    # create folder, if it is not there
    #Path(metadata_emodulus_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    # if config['mode'] == 'cheap' or config['mode'] == 'single':
    #    list_raw_data_emodulus_directories = [ Path(raw_data_emodulus_directory, single_example_name) ]
    # else: # go through all files
    #list_raw_data_emodulus_directories = os.scandir(
        #raw_data_emodulus_directory)

    #for f in list_raw_data_emodulus_directories:
        #if f.is_dir():
            #raw_data_path = Path(f)
            #raw_data_file = Path(f, 'specimen.dat')
            #specimen_metadata_file = Path(
                #metadata_emodulus_directory, f.name + ' _Specimen' + '.json')
            #metadata_file = Path(
                #metadata_emodulus_directory, f.name + '.json')
            #yield {
                #'name': f.name,
                #'actions': [(emodul_metadata, [raw_data_path, metadata_file, specimen_metadata_file])],
                #'file_dep': [raw_data_file],
                #'targets': [metadata_file, specimen_metadata_file],
                #'clean': [clean_targets]
            #}

# extract standardized processed data for Young' modulus tests


#def task_extract_processed_data_emodul():
    # create folder, if it is not there
    #Path(processed_data_emodulus_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    #if config['mode'] == 'cheap' or config['mode'] == 'single':
        #list_raw_data_emodulus_directories = [
            #Path(raw_data_emodulus_directory, single_example_name)]
    #else:  # go through all files
        #list_raw_data_emodulus_directories = os.scandir(
            #raw_data_emodulus_directory)

    #for f in list_raw_data_emodulus_directories:
        #if f.is_dir():
            #raw_data_file = Path(f, 'specimen.dat')
            # the name of the csv file is the file name of the raw data
            # is processed_data_directory + directory_raw_data.csv
            #csv_data_file = Path(
                #processed_data_emodulus_directory, f.name + '.csv')

            #yield {
                #'name': f.name,
                #'actions': [(processed_data_from_rawdata, [f, csv_data_file])],
                #'file_dep': [raw_data_file],
                #'targets': [csv_data_file],
                #'clean': [clean_targets]
            #}

#generate knowledgeGraphs
#@create_after(executed='extract_metadata_emodul')
#def task_export_Mapping():
# create folder, if it is not there
    #Path(knowledge_graphs_output_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    #if config['mode'] == 'cheap' or config['mode'] == 'single':
        #list_metadata_json_files = [Path(metadata_mixture_directory, single_mix_name.split(".xls")[0] + '.json')]
    #else:  # go through all files
        # list of all meta data files....
        #list_metadata_json_files = os.scandir(metadata_mixture_directory)

    # check directory, if

    #for f in list_metadata_json_files:
        #if f.is_file():
            # path to metadata json
            #metadata_file_path = Path(f)
            #name_of_ttl = f.name.replace('.json', '.ttl')
            # path to input KG template
            #KGtemplatePath = Path(knowledge_graphs_Template_directory, single_KG_template)
            # path to output file KG
            #knowledge_graph_file = Path(knowledge_graphs_output_directory, name_of_ttl)

            #yield{
                #'name': f.name,
                #'actions': [(mapping, [KGtemplatePath, metadata_file_path, knowledge_graph_file])],
                #'file_dep': [metadata_file_path],
                #'targets': [knowledge_graph_file],
                #'clean': [clean_targets]
            #}
