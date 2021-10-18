import os
from pathlib import Path
import git

baseDir = Path(__file__).resolve().parents[0]

DOIT_CONFIG = {'verbosity': 2}

def task_create_data_folder():
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
            'actions': ['mkdir E-modul-processed-data/rawdata']
        }
        yield {
            'basename': 'create processeddata folder',
            'actions': ['mkdir E-modul-processed-data/processeddata']
        }

def chowlk_convert_drawio_to_xml():

    def clone_chowlk():
        if os.path.exists(os.path.join(baseDir,'bamChowlk')):
            print('chowlk is already existed')
        else:
            git.Git(baseDir).clone('https://github.com/firmao/bamChowlk')

    yield {
            'actions': [clone_chowlk]
        }