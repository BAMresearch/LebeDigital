from import_gui import main
from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import emodul_metadata

# Script that takes raw data paths as input via GUI and has as output different knowledge graph ttl-files

# start GUI to retrieve paths for raw data as dicitionary, one file for each mixture and for e-module
paths = main()

# run metadata extraction script for mixDesign
mix_metadata(paths['Mixture'], '../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/Workflow2014_08_05 Rezeptur_MI.json')

# run metadata extraction script for e-module (simultaneously creating KG for specimen)
emodul_metadata(paths['Module'], '../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/WorkflowtestMetaData.json','../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/WorkflowtestSpecimenData.json')