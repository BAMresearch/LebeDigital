import subprocess
from time import sleep
from turtle import pos
import requests
import os
import sys
from pathlib import Path
from git import Repo

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
ontologyPath = os.path.join(baseDir2,'ConcreteOntology')
metadataPath = os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv')
graphPath = os.path.join(baseDir0,'E-modul-processed-data/EM_Graph.ttl')
#importedOntologiesPath = os.path.join(baseDir2,'GraphCreation/Ontologies')
processedDataPath = os.path.join(baseDir0,'E-modul-processed-data')

# clone and run rdfconverter docker
try:
    Repo.clone_from('https://github.com/Mat-O-Lab/RDFConverter.git', os.path.join(baseDir1, 'RDFConverter'))
except:
    pass

subprocess.run(['docker-compose', '-f', os.path.join(baseDir1, 'RDFConverter/docker-compose.yml'), 'up', '-d', 'app'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

passed = [] #list of filenames of files that passed validation
failed = [] #list of tuples (filenames, reports) of files that failed validaton

# read shacl shape
with open(os.path.join(baseDir1, 'shape.ttl'), 'r') as f:
    shacl_shape = f.read()

sleep(3)
s = requests.session()

# iterate over all files ending in .ttl
for rdf_file in (os.path.join(processedDataPath, f) for f in os.listdir(processedDataPath) if os.fsdecode(f).endswith('.ttl')):
    with open(rdf_file, 'r') as f:
        post_data = {'rdf': f.read(), 'shaclShape': shacl_shape}
    # post files and obtain result of validation
    r = s.post('http://localhost:8080/rdfconv/validate', files=post_data)
    if r.text.startswith('VALID'):
        passed.append(rdf_file)
    else:
        failed.append((rdf_file, r.text))

# close session object and shutdown docker container
s.close()
subprocess.run(['docker-compose', '-f', os.path.join(baseDir1, 'RDFConverter/docker-compose.yml'), 'down'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

# print passed files
print('files that passed shacl validation:')
for t in enumerate(passed, start=1):
    print(f'{t[0]}. {t[1]}')

# print failed files
print('files that failed shacl validation:')
for t in enumerate(failed, start=1):
    print(f'{t[0]}. {t[1][0]}\n{t[1][1]}')

