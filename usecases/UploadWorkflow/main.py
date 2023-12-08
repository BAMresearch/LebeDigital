import sys
import json
from pathlib import Path

# Fügen Sie das Root-Verzeichnis zum Python-Suchpfad hinzu
root_path = Path(__file__).parent.parent.parent
sys.path.append(str(root_path))

from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import emodul_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_generate_processed_data import processed_data_from_rawdata
from lebedigital.raw_data_processing.Compressive_strength.ComSt_generate_processed_data import processed_rawdata
from lebedigital.raw_data_processing.Compressive_strength.ComSt_metadata_extraction import ComSt_metadata
from lebedigital.mapping.mappingscript import mapping
#from fuseki_upload import upload_ttl_file


# change if necessary
orig_path = 'usecases/UploadWorkflow/files/'
export_path = 'usecases/UploadWorkflow/export/'

# Script that takes raw data paths as input via GUI and has as output different knowledge graph ttl-files
# define output paths for the metadata and paths to Knowledge Graph templates
metadataPaths = {'Mixture_Metadata': 'usecases/UploadWorkflow/files/json',
            'Module_Metadata': 'usecases/UploadWorkflow/files/json',
            'Specimen_EMMetadata': 'usecases/UploadWorkflow/files/json'
                 }

ProcessedDataPath = {'ProcessedData': 'usecases/MinimumWorkingExample/emodul/processed_data',
                     'Processed_Data': 'usecases/MinimumWorkingExample/Druckfestigkeit/processeddata'}

MappedKGtemplate_mix = {'MixMapped': 'usecases/UploadWorkflow/export/Mix.ttl'}
MappedKGtemplate_EModule = {'EModuleMapped': 'usecases/UploadWorkflow/export/EModule.ttl'}
MappedKGtemplate_Specimen = {'SpecimenMapped': 'usecases/UploadWorkflow/export/Specimen.ttl'}

# define input paths to Knowledge Graph templates
KGtemplatePaths = {'Mixture_KG': 'lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl',
                   'Module_KG': 'lebedigital/ConcreteOntology/EModuleOntology_KG_Template.ttl',
                   'Specimen_KG': 'lebedigital/ConcreteOntology/Specimen_KG_Template.ttl',
                   'ComSt_KG': 'lebedigital/ConcreteOntology/CompressiveStrength_KG_Template.ttl'}

rawPaths = {'Mixture-Design': 'usecases/UploadWorkflow/files/mix.xls',
            'E-Module': 'usecases/UploadWorkflow/files/emodule.dat'
}

fuseki_url = "https://mechanics:7fewqg2nDWtZb5HC@fuseki.matolab.org"
dataset_name = "LeBeDigital"

def mix_extraction(institute):
    # if the institute is BAM
    if institute == 'BAM':
        # run metadata extraction script for mixDesign
        path_to_json = mix_metadata(rawPaths['Mixture-Design'], metadataPaths['Mixture_Metadata'])
        # run the mapping script for each extracted metadata to generate three KGs
        mapping(KGtemplatePaths['Mixture_KG'], path_to_json, MappedKGtemplate_mix['MixMapped'])
    else:
        raise ValueError("Error: institute not found")


def emodule_extraction(institute):
    # if the institute is BAM
    if institute == 'BAM':
        # run metadata extraction script for e-module (simultaneously creating KG for specimen)
        processed_data_from_rawdata(rawPaths['E-Module'], ProcessedDataPath['ProcessedData'])
        emodul_metadata(rawPaths['E-Module'], metadataPaths['Module_Metadata'], metadataPaths['Specimen_EMMetadata'],
                        path_to_json)
        mapping(KGtemplatePaths['Specimen_KG'], metadataPaths['Specimen_EMMetadata'],
                MappedKGtemplate_Specimen['SpecimenMapped'])
        mapping(KGtemplatePaths['Module_KG'], metadataPaths['Module_Metadata'],
                MappedKGtemplate_EModule['EModuleMapped'])
    else:
        raise ValueError("Error: institute not found")


def compressive_extraction(institute):
    print("test")


# Open institute.txt and extract institute name
with open(f'{orig_path}institute.txt', 'r') as file:
    institute = file.read()

    # Check if file is empty
    if not institute:
        raise ValueError("Error: institute.txt is empty.")

# Check what files were added in the directory
# file names
file_mix = f'{orig_path}mix.xls'
file_emodule = f'{orig_path}emodule.dat'
file_emodule_xml = f'{orig_path}emodule.xml'
file_compressive = f'{orig_path}compressive.dat'

# Check if files exist
if (Path.cwd() / file_mix).is_file():
    # run function
    mix_extraction(institute)

if (Path.cwd() / file_emodule).is_file() or (Path.cwd() / file_emodule_xml).is_file():
    # run function
    emodule_extraction(institute)

if (Path.cwd() / file_compressive).is_file():
    # run function
    compressive_extraction(institute)
