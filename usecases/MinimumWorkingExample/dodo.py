import os
from pathlib import Path
from doit import create_after
from doit.task import clean_targets

from lebedigital.raw_data_processing.youngs_modulus_data \
    .emodul_metadata_extraction import emodul_metadata

from lebedigital.mapping.emodul_mapping import generate_knowledge_graph

from lebedigital.raw_data_processing.youngs_modulus_data \
    .emodul_generate_processed_data import processed_data_from_rawdata

from doit import get_var

from lebedigital import validation
from lebedigital.validation import SCHEMA

# set a variable to define a cheap or full run
# the default "doit" is set to "doit mode=cheap"
# any other mode value runs the expensive version i.e. "doit mode=full"
config = {"mode": get_var('mode', 'cheap')}

DOIT_CONFIG = {'verbosity': 2}

#parent directory of the minimum working example
ParentDir = os.path.dirname(Path(__file__))

# defining paths
emodul_output_directory = Path(ParentDir, 'emodul')  # folder with metadata yaml files
raw_data_emodulus_directory = Path(ParentDir, 'Data', 'E-modul')  # folder with folders of raw data files
metadata_emodulus_directory = Path(emodul_output_directory, 'metadata_yaml_files')  # folder with metadata yaml files
processed_data_emodulus_directory = Path(emodul_output_directory, 'processed_data')  # folder with csv data files
knowledge_graphs_directory = Path(emodul_output_directory, 'knowledge_graphs')  # folder with KG ttl files

# when "cheap option" is run, only this souce of raw data is processed
cheap_example_name = 'Wolf 8.2 Probe 1'

# create folder, if it is not there
Path(emodul_output_directory).mkdir(parents=True, exist_ok=True)

#extract standardized meta data for Young' modulus tests
def task_extract_metadata_emodul():
    # create folder, if it is not there
    Path(metadata_emodulus_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap':
        list_raw_data_emodulus_directories = [ Path(raw_data_emodulus_directory, cheap_example_name) ]
    else: # go through all files
        list_raw_data_emodulus_directories = os.scandir(raw_data_emodulus_directory)

    for f in list_raw_data_emodulus_directories:
        if f.is_dir():
            raw_data_path = Path(f)
            raw_data_file = Path(f, 'specimen.dat')
            yaml_metadata_file = Path(metadata_emodulus_directory, f.name + '.yaml')
            yield {
                'name': yaml_metadata_file,
                'actions': [(emodul_metadata, [raw_data_path, yaml_metadata_file])],
                'file_dep': [raw_data_file],
                'targets': [yaml_metadata_file],
                'clean': [clean_targets]
            }

#extract standardized processed data for Young' modulus tests
def task_extract_processed_data_emodul():
    # create folder, if it is not there
    Path(processed_data_emodulus_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap':
        list_raw_data_emodulus_directories = [ Path(raw_data_emodulus_directory, cheap_example_name) ]
    else: # go through all files
        list_raw_data_emodulus_directories = os.scandir(raw_data_emodulus_directory)

    for f in list_raw_data_emodulus_directories:
        if f.is_dir():
            raw_data_file = Path(f, 'specimen.dat')
            #the name of the csv file is the file name of the raw data
            # is processed_data_directory + directory_raw_data.csv
            csv_data_file = Path(processed_data_emodulus_directory, f.name + '.csv')

            yield {
                'name': csv_data_file,
                'actions': [(processed_data_from_rawdata, [f, csv_data_file])],
                'file_dep': [raw_data_file],
                'targets': [csv_data_file],
                'clean': [clean_targets]
            }


#generate knowledgeGraphs
@create_after(executed='extract_metadata_emodul')
def task_export_knowledgeGraph_emodul():
    # create folder, if it is not there
    Path(knowledge_graphs_directory).mkdir(parents=True, exist_ok=True)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap':
        list_metadata_yaml_files = [ Path(metadata_emodulus_directory, cheap_example_name + '.yaml') ]
    else: # go through all files
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
            processed_data_file_path = Path(processed_data_emodulus_directory, name_of_cvs)
            # path to output file KG
            knowledge_graph_file = Path(knowledge_graphs_directory, name_of_ttl)

            yield{
                'name': knowledge_graph_file,
                'actions': [(generate_knowledge_graph, [metadata_file_path,
                                                    knowledge_graph_file])],
                'file_dep': [metadata_file_path, processed_data_file_path],
                'targets': [knowledge_graph_file],
                'clean': [clean_targets]
            }

#validate
@create_after(executed='expord_knowledgeGraph_emodul')
def task_validate_graph():
    
    graphs = os.scandir(knowledge_graphs_directory)

    s = validation.read_graph_from_file(Path(emodul_output_directory, 'validation_test', 'shape.ttl'))
    
    for f in graphs:
        if f.is_file() and Path(f).suffix == '.ttl':
            # do some validation
            g = validation.read_graph_from_file(g)
            res = validation.test_graph(g, s)

            assert not validation.violates_shape(res, SCHEMA.SpecimenDiameterShape)
            assert not validation.violates_shape(res, SCHEMA.SpecimenShape)
            assert validation.violates_shape(res, SCHEMA.InformationBearingEntityShape)

            out = open(Path(emodul_output_directory, 'validation_result.txt'), 'wx')
            out.write(f'{f.name}:')

            for shape in [
                SCHEMA.SpecimenDiameterShape, 
                SCHEMA.SpecimenShape, 
                SCHEMA.InformationBearingEntityShape
                ]:
                out.write(repr(shape) + ('failed' if validation.violates_shape(res, shape) else 'passed'))
            