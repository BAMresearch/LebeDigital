import os
from pathlib import Path

BASEDIR = Path(__file__).resolve().parents[0]
BASEDIR1 = Path(__file__).resolve().parents[0]
EMODUL_FOLDER = os.path.join(BASEDIR, 'emodul')
EMODUL_RAWDATA_FOLDER = os.path.join(EMODUL_FOLDER,'rawdata')
EMODUL_PROCESSEDDATA_FOLDER = os.path.join(EMODUL_FOLDER,'processeddata')
EMODUL_YAML_METADATA_FOLDER = os.path.join(EMODUL_FOLDER,'metadata_yaml_files')

EMODUL_METADATA_EXTRACTION_SCRIPT = os.path.join(os.path.join(os.path.join(BASEDIR1, 'src'),'metadata_extraction'),'emodul_metadata_extraction.py')

DOIT_CONFIG = {'verbosity': 2}

def task_installation():
    yield {
        'basename': 'install python packages',
        'actions': ['pip install -r requirements.txt'],
        'file_dep': ['requirements.txt']
    }

def task_emodul():
    if os.path.exists(os.path.join(BASEDIR,'emodul')):
        yield {
            'basename': 'checking folder',
            'actions': ['echo folder existed']
        }
    else:

        yield {
            'basename': 'create E-modul-processed-data folder',
            'actions': ['mkdir %s ' % EMODUL_FOLDER]
        }
        yield {
            'basename': 'create rawdata folder',
            'actions': ['mkdir %s ' % EMODUL_RAWDATA_FOLDER]
        }
        yield {
            'basename': 'create processeddata folder',
            'actions': ['mkdir %s ' % EMODUL_PROCESSEDDATA_FOLDER]
        }
        yield {
            'basename': 'create folder contains metadata yaml files',
            'actions': ['mkdir %s ' % EMODUL_YAML_METADATA_FOLDER]
        }
    
    yield {
        'basename': 'generate emodul processed data',
        'actions': ['python %s' % EMODUL_METADATA_EXTRACTION_SCRIPT],
        # 'file_dep': ['knowledgeGraph/emodul/emodul_generate_processed_data.py']
    }
#     yield {
#         'basename': 'extract emodul metadata',
#         'actions': ['python knowledgeGraph/emodul/emodul_metadata_extraction.py'],
#         'targets': ['knowledgeGraph/emodul/E-modul-processed-data/emodul_metadata.csv'],
#         'file_dep': ['knowledgeGraph/emodul/emodul_metadata_extraction.py']
#     }
#     yield {
#         'basename': 'calculate emodul',
#         'actions': ['python knowledgeGraph/emodul/emodul_calculation.py'],
#         'file_dep': ['knowledgeGraph/emodul/emodul_calculation.py']
#     }
#     yield {
#         'basename': 'map emodul ontology and metadata',
#         'actions': ['python knowledgeGraph/emodul/emodul_mapping.py'],
#         'targets': ['knowledgeGraph/emodul/E-modul-processed-data/EM_Graph.ttl'],
#         'file_dep': ['knowledgeGraph/emodul/emodul_mapping.py']
#     }
#     # yield {
#     #     'basename': 'validate rdf files against shacl shape',
#     #     'actions': ['python knowledgeGraph/emodul/emodul_validation.py']
#     # }
#     yield {
#         'basename': 'run emodul query script',
#         'actions': ['python knowledgeGraph/emodul/emodul_query.py'],
#         'file_dep': ['knowledgeGraph/emodul/emodul_query.py']
#     }
#     yield {
#         'basename': 'run emodul test query',
#         'actions': ['python knowledgeGraph/emodul/emodul_test.py'],
#         'file_dep': ['knowledgeGraph/emodul/emodul_test.py']
#     }

