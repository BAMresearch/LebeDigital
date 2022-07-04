import os
from pathlib import Path
import sys
from doit import create_after
from doit.tools import create_folder

from lebedigital.raw_data_processing.metadata_extraction \
    .emodul_metadata_extraction import eModul_metadata

from lebedigital.raw_data_processing.processed_data_generation \
    .emodul_generate_processed_data import processed_data_from_rawdata


DOIT_CONFIG = {'verbosity': 2}
#parent directory of the minimum working example
ParentDir = os.path.dirname(Path(__file__))

#raw data files is a subdirectory
raw_data_emodulus_directory = Path(ParentDir, 'Data', 'E-modul')
#meta data is in a different subdirectory
meta_data_emodulus_directory = Path(ParentDir, 'emodul',
                                    'meta_data_yaml_files')
#processed data is in a different subdirectory
processed_data_emodulus_directory = Path(ParentDir, 'emodul',
                                    'processed_data')

#extract standardized meta data for Young' modulus tests
def task_extract_metadata_emodul():
    for f in os.scandir(raw_data_emodulus_directory):
        if f.is_dir():
            raw_data_file = Path(f, 'specimen.dat')
            yaml_meta_data_file = Path(meta_data_emodulus_directory, f.name + '.yaml')
            yield {
                'name': yaml_meta_data_file,
                'actions': [(eModul_metadata, [raw_data_file,
                                               meta_data_emodulus_directory,
                                               f.name + '.yaml'])],
                'file_dep': [raw_data_file],
                'targets': [yaml_meta_data_file],
            }


#extract standardized processed data for Young' modulus tests
def task_extract_processed_data_emodul():
    for f in os.scandir(raw_data_emodulus_directory):
        if f.is_dir():
            raw_data_file = Path(f, 'specimen.dat')
            #the name of the csv file is the file name of the raw data
            # is processed_data_directory + directory_raw_data.csv
            csv_data_file = Path(processed_data_emodulus_directory,
                                 f.name + '.csv')
            yield {
                'name': csv_data_file,
                'actions': [(processed_data_from_rawdata, [raw_data_file,
                                               processed_data_emodulus_directory,
                                               f.name + '.csv'])],
                'file_dep': [raw_data_file],
                'targets': [csv_data_file],
            }