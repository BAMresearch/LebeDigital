import os
from pathlib import Path

baseDir = Path(__file__).resolve().parents[0]
emodulRawdataFolder = os.path.join(os.path.join(os.path.join(baseDir,'emodul'),'E-modul-processed-data'),'rawdata')
emodulProcesseddataFolder = os.path.join(os.path.join(os.path.join(baseDir,'emodul'),'E-modul-processed-data'),'processeddata')

compressionRawdataFolder = os.path.join(os.path.join(os.path.join(baseDir,'compression'),'compression-processed-data'),'rawdata')
compressionProcesseddataFolder = os.path.join(os.path.join(os.path.join(baseDir,'compression'),'compression-processed-data'),'processeddata')

DOIT_CONFIG = {'verbosity': 2}

def task_installation():
    yield {
        'basename': 'install python packages',
        'actions': ['pip install -r requirements.txt']
    }

def task_emodul():
    if os.path.exists(os.path.join(os.path.join(baseDir,'emodul'),'E-modul-processed-data')):
        yield {
            'basename': 'checking folder',
            'actions': ['echo folder existed']
        }
    else:

        yield {
            'basename': 'create E-modul-processed-data folder',
            'actions': ['mkdir emodul/E-modul-processed-data']
        }
        yield {
            'basename': 'create rawdata folder',
            'actions': ['mkdir %s ' % emodulRawdataFolder]
        }
        yield {
            'basename': 'create processeddata folder',
            'actions': ['mkdir %s ' % emodulProcesseddataFolder]
        }
    
    yield {
        'basename': 'generate emodul processed data',
        'actions': ['python emodul/emodul_generate_processed_data.py']
    }
    yield {
        'basename': 'extract emodul metadata',
        'actions': ['python emodul/emodul_metadata_extraction.py']
    }
    yield {
        'basename': 'map emodul ontology and metadata',
        'actions': ['python emodul/emodul_mapping.py']
    }
    yield {
        'basename': 'run emodul query script',
        'actions': ['python emodul/emodul_query.py']
    }
    yield {
        'basename': 'run emodul test query',
        'actions': ['python emodul/emodul_test.py']
    }

def task_compression():
    if os.path.exists(os.path.join(os.path.join(baseDir,'compression'),'compression-processed-data')):
        yield {
            'basename': 'checking compression processed data folder',
            'actions': ['echo folder existed']
        }
    else:

        yield {
            'basename': 'create compression-processed-data folder',
            'actions': ['mkdir compression/compression-processed-data']
        }
        yield {
            'basename': 'create compression rawdata folder',
            'actions': ['mkdir %s ' % compressionRawdataFolder]
        }
        yield {
            'basename': 'create compression processeddata folder',
            'actions': ['mkdir %s ' % compressionProcesseddataFolder]
        }
    yield {
        'basename': 'generate compression processed data',
        'actions': ['python compression/compression_generate_processed_data.py']
    }
    yield {
        'basename': 'extract compression metadata',
        'actions': ['python compression/compression_metadata_extraction.py']
    }
    yield {
        'basename': 'map compression ontology and metadata',
        'actions': ['python compression/compression_mapping.py']
    }
    yield {
        'basename': 'run compression query script',
        'actions': ['python compression/compression_query.py']
    }
    yield {
        'basename': 'run compression test query',
        'actions': ['python compression/compression_test.py']
    }