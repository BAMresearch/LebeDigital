from import_gui import main
from Lebedigital.lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from Lebedigital.lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import emodul_metadata
from Lebedigital.lebedigital.mapping.mappingscript import mapping

# Script that takes raw data paths as input via GUI and has as output different knowledge graph ttl-files
# define output paths for the metadata and paths to Knowledge Graph templates
metadataPaths = {'Mixture_Metadata': '../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/',
            'Module_Metadata': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/testMetaData.json',
            'Specimen_Metadata': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/testSpecimenData.json'}


MappedKGtemplate_mix = {'MixMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMixMapped.ttl'}
MappedKGtemplate_EModule = {'EModuleMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMapped.ttl'}
MappedKGtemplate_Specimen = {'SpecimenMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testSpecimenMapped.ttl'}

# define input paths to Knowledge Graph templates
KGtemplatePaths = {'Mixture_KG': '../../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl',
                   'Module_KG': '../../../lebedigital/ConcreteOntology/EModuleOntology_KG_Template.ttl',
                   'Specimen_KG': '../../../lebedigital/ConcreteOntology/Specimen_KG_Template.ttl'}

# start GUI to retrieve paths for raw data as dictionary, one file for each mixture and for e-module
rawPaths = main()

# run metadata extraction script for mixDesign
path_to_json = mix_metadata(rawPaths['Mixture Design'], metadataPaths['Mixture_Metadata'])

# run metadata extraction script for e-module (simultaneously creating KG for specimen)
emodul_metadata(rawPaths['E-Module'], metadataPaths['Module_Metadata'], metadataPaths['Specimen_Metadata'], path_to_json)

# run the mapping script for each extracted metadata to generate three KGs
#mapping(KGtemplatePaths['Mixture_KG'], metadataPaths['Mixture_Metadata'], MappedKGtemplate_mix['MixMapped'])
mapping(KGtemplatePaths['Mixture_KG'], path_to_json, MappedKGtemplate_mix['MixMapped'])
mapping(KGtemplatePaths['Specimen_KG'], metadataPaths['Specimen_Metadata'], MappedKGtemplate_Specimen['SpecimenMapped'])
mapping(KGtemplatePaths['Module_KG'], metadataPaths['Module_Metadata'], MappedKGtemplate_EModule['EModuleMapped'])



