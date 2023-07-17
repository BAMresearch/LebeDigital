from import_gui import main
from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import emodul_metadata
from lebedigital.mapping.mappingscript import mapping

# Script that takes raw data paths as input via GUI and has as output different knowledge graph ttl-files

# define output paths for the metadata and paths to Knowledge Graph templates
metadataPaths = {'Mixture': '../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/Workflow2014_08_05 Rezeptur_MI.json',
            'Module:' : '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/WorkflowtestMetaData.json',
            'Specimen': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/WorkflowtestSpecimenData.json'}

KGtemplatePaths = {'Mixture' : '../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl',
                   'Module' : '../../lebedigital/ConcreteOntology/EModuleOntology_KG_Template.ttl',
                   'Specimen' : '../../lebedigital/ConcreteOntology/Specimen_KG_Template.ttl'}

# start GUI to retrieve paths for raw data as dictionary, one file for each mixture and for e-module
rawPaths = main()

# run metadata extraction script for mixDesign
mix_metadata(rawPaths['Mixture'], metadataPaths['Mixture'])

# run metadata extraction script for e-module (simultaneously creating KG for specimen)
emodul_metadata(rawPaths['Module'], metadataPaths['Module'], metadataPaths['Specimen'])

# run the mapping script for each extracted metadata
mapping(metadataPaths['Mixture'], KGtemplatePaths['Mixture'])
mapping(metadataPaths['Module'], KGtemplatePaths['Module'])
mapping(metadataPaths['Specimen'], KGtemplatePaths['Specimen'])