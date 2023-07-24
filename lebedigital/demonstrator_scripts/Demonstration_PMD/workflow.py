from import_gui import main
from Lebedigital.lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from Lebedigital.lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import emodul_metadata
from Lebedigital.lebedigital.mapping.mappingscript import mapping

# Script that takes raw data paths as input via GUI and has as output different knowledge graph ttl-files

# define output paths for the metadata and paths to Knowledge Graph templates
metadataPaths = {'Mixture_Metadata': '../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/Workflow2014_08_05 Rezeptur_MI.json',
            'Module_Metadata': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/WorkflowtestMetaData.json',
            'Specimen_Metadata': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/WorkflowtestSpecimenData.json'}


MappedKGtemplate_mix = {'MixMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMixMapped.ttl'}
MappedKGtemplate_EModule = {'EModuleMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMapped.ttl'}
MappedKGtemplate_Specimen = {'SpecimenMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMapped.ttl'}

# define input paths to Knowledge Graph templates
KGtemplatePaths = {'Mixture_KG': '../../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl',
                   'Module_KG': '../../../lebedigital/ConcreteOntology/EModuleOntology_KG_Template.ttl',
                   'Specimen_KG': '../../../lebedigital/ConcreteOntology/Specimen_KG_Template.ttl'}

# start GUI to retrieve paths for raw data as dictionary, one file for each mixture and for e-module
rawPaths = main()

# run metadata extraction script for mixDesign
mix_metadata(rawPaths['Mixture Design'], metadataPaths['Mixture_Metadata'])

# run metadata extraction script for e-module (simultaneously creating KG for specimen)
emodul_metadata(rawPaths['E-Module'], metadataPaths['Module_Metadata'], metadataPaths['Specimen_Metadata'])

# run the mapping script for each extracted metadata to generate three KGs
mapping(KGtemplatePaths['Mixture_KG'], metadataPaths['Mixture_Metadata'], MappedKGtemplate_mix['MixMapped'])
mapping(metadataPaths['Module_Metadata'], KGtemplatePaths['Module_KG'], MappedKGtemplate_EModule['EModuleMapped'])
mapping(metadataPaths['Specimen_Metadata'], KGtemplatePaths['Specimen_KG'], MappedKGtemplate_Specimen['SpecimenMapped'])