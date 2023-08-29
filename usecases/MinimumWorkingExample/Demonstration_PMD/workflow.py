from import_gui import main
from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import emodul_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_generate_processed_data import processed_data_from_rawdata
from lebedigital.mapping.mappingscript import mapping
from Merge_KGs import merge_graphs
from fuseki_upload import upload_ttl_file

# Script that takes raw data paths as input via GUI and has as output different knowledge graph ttl-files
# define output paths for the metadata and paths to Knowledge Graph templates
metadataPaths = {'Mixture_Metadata': '../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/',
            'Module_Metadata': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/testMetaData.json',
            'Specimen_Metadata': '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/testSpecimenData.json'}

ProcessedDataPath = {'ProcessedData': '../../../usecases/MinimumWorkingExample/emodul/processed_data'}

MappedKGtemplate_mix = {'MixMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMixMapped.ttl'}
MappedKGtemplate_EModule = {'EModuleMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testMapped.ttl'}
MappedKGtemplate_Specimen = {'SpecimenMapped': '../../../usecases/MinimumWorkingExample/Mapping_Example/testSpecimenMapped.ttl'}



# define input paths to Knowledge Graph templates
KGtemplatePaths = {'Mixture_KG': '../../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl',
                   'Module_KG': '../../../lebedigital/ConcreteOntology/EModuleOntology_KG_Template.ttl',
                   'Specimen_KG': '../../../lebedigital/ConcreteOntology/Specimen_KG_Template.ttl'}

fuseki_url = "https://mechanics:7fewqg2nDWtZb5HC@fuseki.matolab.org"
dataset_name = "LeBeDigital"

# start GUI to retrieve paths for raw data as dictionary, one file for each mixture and for e-module
rawPaths = main()

# run metadata extraction script for mixDesign
path_to_json = mix_metadata(rawPaths['Mixture Design'], metadataPaths['Mixture_Metadata'])

# run metadata extraction script for e-module (simultaneously creating KG for specimen)
emodul_metadata(rawPaths['E-Module'], metadataPaths['Module_Metadata'], metadataPaths['Specimen_Metadata'], path_to_json)
processed_data_from_rawdata(rawPaths['E-Module'], ProcessedDataPath['ProcessedData'])

# run the mapping script for each extracted metadata to generate three KGs
mapping(KGtemplatePaths['Mixture_KG'], path_to_json, MappedKGtemplate_mix['MixMapped'])
mapping(KGtemplatePaths['Specimen_KG'], metadataPaths['Specimen_Metadata'], MappedKGtemplate_Specimen['SpecimenMapped'])
mapping(KGtemplatePaths['Module_KG'], metadataPaths['Module_Metadata'], MappedKGtemplate_EModule['EModuleMapped'])

# Merge ttls and upload it on fuseki
merged_graph = merge_graphs(MappedKGtemplate_mix['MixMapped'], MappedKGtemplate_Specimen['SpecimenMapped'], MappedKGtemplate_EModule['EModuleMapped'])
upload_ttl_file(fuseki_url, dataset_name, merged_graph)
