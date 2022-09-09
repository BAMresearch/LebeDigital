import os
from pathlib import Path
from doit import create_after
from doit.task import clean_targets

from lebedigital.raw_data_processing.metadata_extraction \
    .emodul_metadata_extraction import emodul_metadata

from lebedigital.mapping.emodul_mapping import generate_knowledge_graph

from lebedigital.raw_data_processing.processed_data_generation \
    .emodul_generate_processed_data import processed_data_from_rawdata
from lebedigital.calibration.calibrationWorkflow import perform_calibration

from doit import get_var
from doit.tools import config_changed

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
calibrated_data_directory = Path(emodul_output_directory, 'calibrated_data')  # folder with calibration output

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
                'clean': [clean_targets]  # what does this do?
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

# perform calibration
def task_perform_calibration():
    # create folder, if it is not there
    Path(calibrated_data_directory).mkdir(parents=True, exist_ok=True)



    # TODO implement "expensive" option on all files!!!
    knowledge_graphs_file = Path(knowledge_graphs_directory, cheap_example_name + '.ttl')

    # current outputs
    owl_file = Path(calibrated_data_directory, 'calibrationWorkflow.owl')
    displ_list = Path(calibrated_data_directory, 'displacement_list_' + cheap_example_name + '.dat')
    force_list = Path(calibrated_data_directory, 'force_list_' + cheap_example_name + '.dat')
    another_file = Path(calibrated_data_directory, 'joint_samples_compression_test_calibration.dat')

    yield {
        'name': f'calibrate {cheap_example_name}',
        'actions': [(perform_calibration,[str(knowledge_graphs_directory),calibrated_data_directory,cheap_example_name,config['mode']])],
        'file_dep': [knowledge_graphs_file],
        'uptodate': [config_changed(config['mode'])],  # the results of the calibration depend on the mode, therfore it must be recomputed when the mode changes
        'targets': [owl_file,displ_list,force_list,another_file],
        'clean': [clean_targets]
    }