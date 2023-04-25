from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import extract_metadata_emodulus
from pathlib import Path

def test_metadata_extraction_emodul():
        
        """Testing the emodule metadata extraction on a single example (2014_08_05 Rezeptur_MI)"""

        # setup paths and directories
        data_dir = 'test_data'
        data_path = Path(__file__).parent / data_dir

        mix_file = 'mix.dat'
        specimen_file = 'specimen.dat'

        target_data = {"CompressionColumn": 0,
                        "ExperimentDate": "01.09.2014",
                        "ExperimentName": "test_data",
                        "ExperimentTime": '13:25:39',
                        "Lab": "BAM",
                        "MeasurementDuration": 320.02344,
                        "MixDataFile": "../LebeDigital/tests/raw_data_processing/Mischungen/2014_08_05 Rezeptur_MI.xlsx",
                        "ProcessedFile": "../usecases/MinimumWorkingExample/emodul/processed_data",
                        "RawDataFile": "../LebeDigital/tests/raw_data_processing/youngs_modulus_data/test_data/specimen.dat",
                        "SpecimenAge": 28.0,
                        "SpecimenDiameter": 98.6,
                        "SpecimenLength": 300.3,
                        "SpecimenName": "BA-Losert E-Modul 28d v. 04.08.14 Probe 4",
                        "SpecimenWeight": 5342.0,
                        "TransducerColumn": [1,2,3],
                        "remark": "Kraftgeregelt 3,9 kN/s",
                        "software_specification": "MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A",
                        "tester_name": "Kh"
                        }

        # run extraction and getting a dictionary with metadata
        test_data = extract_metadata_emodulus(data_path, specimen_file, mix_file)

        assert test_data == target_data
  
