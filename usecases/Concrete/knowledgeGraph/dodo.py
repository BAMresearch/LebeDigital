import os
from pathlib import Path

baseDir = Path(__file__).resolve().parents[0]
rawdataFolder = os.path.join(os.path.join(baseDir,'E-modul-processed-data'),'rawdata')
processeddataFolder = os.path.join(os.path.join(baseDir,'E-modul-processed-data'),'processeddata')

DOIT_CONFIG = {'verbosity': 2}

def task_emodul():
    if os.path.exists(os.path.join(baseDir,'E-modul-processed-data')):
        yield {
            'basename': 'checking folder',
            'actions': ['echo folder existed']
        }
    else:

        yield {
            'basename': 'create E-modul-processed-data folder',
            'actions': ['mkdir E-modul-processed-data']
        }
        yield {
            'basename': 'create rawdata folder',
            'actions': ['mkdir %s ' % rawdataFolder]
        }
        yield {
            'basename': 'create processeddata folder',
            'actions': ['mkdir %s ' % processeddataFolder]
        }
    yield {
        'basename': 'install python packages',
        'actions': ['pip install -r requirements.txt']
    }

    yield {
        'basename': 'generate processed data',
        'actions': ['python generate_Emodul_processed_data.py']
    }
    yield {
        'basename': 'extract metadata',
        'actions': ['python metadata_extraction.py']
    }
    yield {
        'basename': 'map ontology and metadata',
        'actions': ['python mapping_ontology_data.py']
    }
    yield {
        'basename': 'run query script',
        'actions': ['python query.py']
    }

    