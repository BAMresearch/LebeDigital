from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import extract_metadata_emodulus
from pathlib import Path

def test_metadata_extraction_emodul():
        """Tesing the metadata extraction on a single example"""

        # setup paths and directories
        data_dir = 'test_data'
        data_path = Path(__file__).parent / data_dir

        mix_file = 'mix.dat'
        specimen_file = 'specimen.dat'

        target_data = {
                         "humanreadableID": "BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4",
                         "ID": "d8ce15fe-bd5c-446b-baff-cb44d970e16a",
                         "MixtureID": "6fea0bac-e04b-458d-b7b5-f8e2712c5efa",
                         "SpecimenMass": 5342.0,
                         "SpecimenMass_Unit": "g",
                         "SpecimenDiameter": 98.6,
                         "SpecimenDiameter_Unit": "mm",
                         "SpecimenLength": 300.3,
                         "SpecimenLength_Unit": "mm"
                      }
                      

        # run extraction and getting a dictionary with metadata
        test_metadata_emodule, test_metadata_specimen = extract_metadata_emodulus(data_path, specimen_file, mix_file)

        assert test_metadata_specimen == target_data
