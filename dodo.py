import os
from scripts.mapping.mappingscript import mapping
from scripts.mapping.mixmapping import mappingmixture
from scripts.upload.upload_script import upload_to_docker
from scripts.upload.json_db import upload_to_db

# bearer token, need to find a way to make it private
DOIT_CONFIG = {
    'docker_bearer': ""
}

# Paths
import_path = 'files/json'
export_path = 'files/ttl'
kg_templates_paths = {'Mixture_KG': 'cpto/MixtureDesign_KG_Template.ttl',
                   'Module_KG': 'cpto/EModuleOntology_KG_Template.ttl',
                   'Specimen_KG': 'cpto/Specimen_KG_Template.ttl',
                   'ComSt_KG': 'cpto/CompressiveStrength_KG_Template.ttl'}
output_files_name = {'Mixture': 'Mixture_Mapped.ttl',
                     'Module': 'Module.ttl',
                     'Specimen': 'Specimen.ttl',
                     'Compressive': 'Compressive.ttl'}

# list of all datasets that were sucessfully mapped, to be uploaded to ontodocker
datasets = []


def task_mix_mapping():
    def mix_mapping():
        dateipfad = f'{import_path}/mix.json'
        # check if the file exists
        if not os.path.isfile(dateipfad):
            print(f"Die Datei {dateipfad} existiert nicht.")
            return

        mappingmixture(dateipfad, f'{export_path}/{output_files_name["Mixture"]}')
        datasets.append(f'{export_path}/{output_files_name["Mixture"]}')
        return
    return {
        'actions': [mix_mapping],
        'verbosity': 2,
    }


def task_em_mapping():
    def em_mapping():
        dateipfad = f'{import_path}/emodule.json'
        dateipfad2 = f'{import_path}/specimen.json'
        # check if the file exists
        if not os.path.isfile(dateipfad) or not os.path.isfile(dateipfad2):
            print(f"Die Datei {dateipfad} oder {dateipfad2} existiert nicht.")
            return

        mapping(kg_templates_paths['Specimen_KG'], dateipfad2, f'{export_path}/{output_files_name["Specimen"]}')
        mapping(kg_templates_paths['Module_KG'], dateipfad, f'{export_path}/{output_files_name["Module"]}')
        datasets.append(f'{export_path}/{output_files_name["Specimen"]}')
        datasets.append(f'{export_path}/{output_files_name["Module"]}')
        return
    return {
        'actions': [em_mapping],
        'verbosity': 2,
    }


def task_upload_to_ontodocker():
    def uploadtoontodocker():
        upload_to_docker(DOIT_CONFIG['docker_bearer'], datasets)
        return
    return {
        'actions': [uploadtoontodocker],
        'verbosity': 2,
    }

def task_clean_up():
    def clean_up():
        rem_files = ['files/ttl/Specimen.ttl', 'files/ttl/Module.ttl', 'files/ttl/Mixture_Mapped.ttl', 'files/json/emodule.json', 'files/json/mix.json', 'files/json/specimen.json']
        # remove ttl files
        for entry in rem_files:
            os.remove(entry)
        print('Clean up complete.')
    return {
        'actions': [clean_up],
        'task_dep': ['json_to_db', 'mix_mapping', 'em_mapping', 'upload_to_ontodocker'],
        'verbosity': 2,
    }