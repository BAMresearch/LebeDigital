from doctest import Example
import os
from pathlib import Path
import sys
from doit import create_after
from doit.tools import create_folder

DOIT_CONFIG = {'verbosity': 2}
PARENTDIR = os.path.dirname(Path(__file__).resolve().parents[0])
sys.path.insert(0, PARENTDIR)
from lebedigital.raw_data_processing.metadata_extraction\
    .emodul_metadata_extraction import eModul_metadata
from lebedigital.raw_data_processing.processed_data_generation.emodul_generate_processed_data import processed_data_from_rawdata
from lebedigital.mapping.emodul_mapping import metadata_ontology_mapping

emodulRawDataFolder = os.path.join(os.path.join(os.path.join(PARENTDIR, 'Example'), 'Data'),'E-modul')
emodulResultFolder = os.path.join(os.path.join(PARENTDIR, 'Example'), 'emodul')
emodulMetadataFolder = os.path.join(os.path.join(os.path.join(PARENTDIR, 'Example'), 'emodul'), 'metadata_yaml_files')
emodulProcessedDataFolder = os.path.join(os.path.join(os.path.join(PARENTDIR, 'Example'), 'emodul'), 'processeddata')
emodulTriplesFolder = os.path.join(os.path.join(os.path.join(PARENTDIR, 'Example'), 'emodul'), 'triples')


def task_emodul():
    def emodul_metadata_extraction(locationOfRawData, fileName,targets):
        eModul_metadata(locationOfRawData, fileName,targets[0])

    listRawData = os.listdir(emodulRawDataFolder)
    for rawData in listRawData:
        yield {
            'actions': [(emodul_metadata_extraction, [os.path.join(emodulRawDataFolder, rawData), 'specimen.dat'])],
            'targets': [os.path.join(emodulMetadataFolder, rawData + '.yaml')],
            'basename': rawData + ' metadata extraction'
        }

    def emodul_processed_data_generation(locationOfRawData,targets):
        processed_data_from_rawdata(locationOfRawData,targets[0])

    listRawData = os.listdir(emodulRawDataFolder)
    for rawData in listRawData:
        yield {
            'actions': [(emodul_processed_data_generation, [os.path.join(emodulRawDataFolder, rawData)])],
            'targets': [os.path.join(emodulProcessedDataFolder, 'processed_' + rawData.replace('.','_').replace(' ','_') + '.csv')],
            'basename': rawData + ' processed data generation'
        }
    
    def emodul_knowledge_graph_generation(locationOfMetadata,targets):
        metadata_ontology_mapping(locationOfMetadata,targets[0])

    listMetadata = os.listdir(emodulMetadataFolder)
    for metadata in listMetadata:
        yield {
            'actions': [(emodul_knowledge_graph_generation, [os.path.join(emodulMetadataFolder, metadata)])],
            'targets': [os.path.join(emodulTriplesFolder, metadata.replace('.yaml', '.ttl'))],
            'basename': metadata.replace('.yaml', '') + ' knowledge graph generation',
        }

