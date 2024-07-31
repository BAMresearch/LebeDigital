import os
import json
from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import mix_metadata
from pathlib import Path


def test_mixdesign_metadata_extraction():
    """Testing the mixture metadata extraction on a single example"""

    # setup paths and directories
    data_dir = 'test_data'
    data_path = Path(__file__).parent / data_dir
    input_file = '20240220_7188_M01.xls'
    input_path = data_path / input_file
    output_file = 'test_'
    output_path = data_path / (output_file + '.json')  # Assuming output is stored as a JSON file

    # make sure that the path exists
    assert os.path.exists(input_path)

    target_data = {
        "RawDataFile": "test_data/20240220_7188_M01.xls",
        "MixingDate": "2024-02-20T12:00:00",
        "Lab": "BAM",
        "ID": "4e04a345-132a-44ed-ae60-6e065c201d2f",
        "humanreadableID": "20240220_7188_M01",
        "Cement1_Content": 375.0,
        "Cement1_Content_Unit": "kg/m^3",
        "Cement1_Density": 3.1,
        "Cement1_Density_Unit": "kg/dm^3",
        "Cement1_Type": "CEM I 42.5 N Rüdersdorf",
        "Water_Content": 168.75,
        "Water_Content_Unit": "kg/m^3",
        "Water_Density": None,
        "Water_Density_Unit": "kg/dm^3",
        "Water_Type": "incl. Wasser aus Fließmittel (73.4 M.-%)",
        "WaterCementRatio": 0.45,
        "Admixture1_Content": 5.625,
        "Admixture1_Content_Unit": "kg/m^3",
        "Admixture1_Density": 1.13,
        "Admixture1_Density_Unit": "kg/dm^3",
        "Admixture1_Type": "MasterRheobuild 1021 Fließmittel",
        "Admixture2_Content": 0.0,
        "Admixture2_Content_Unit": "kg/m^3",
        "Admixture2_Density": 1.01,
        "Admixture2_Density_Unit": "kg/dm^3",
        "Admixture2_Type": "LP-Mittel",
        "Aggregate1_Content": 1830.0,
        "Aggregate1_Content_Unit": "kg/m^3",
        "Aggregate1_Size": None,
        "Aggregate1_Size_Unit": None,
        "Aggregate1_Density": None,
        "Aggregate1_Density_Unit": "kg/dm^3",
        "Addition1_Content": 0.0,
        "Addition1_Content_Unit": "kg/m^3",
        "Addition1_Density": 2.33,
        "Addition1_Density_Unit": "kg/dm^3",
        "Addition1_Type": "gesamt"
    }

    # Keys to ignore during comparison
    keys_to_ignore = ["RawDataFile", "ID"]

    # Run extraction and save to output JSON file
    test_metadata = mix_metadata(input_path, output_file)

    # Load the JSON file
    with open(test_metadata + '.json', 'r') as file:
        test_data = json.load(file)

    # Debugging print statements
    print("Type of test_data:", type(test_data))
    print("Content of test_data:", test_data)

    # Create dictionaries without the keys to ignore
    target_data_without_keys = {key: value for key, value in target_data.items() if key not in keys_to_ignore}
    test_data_without_keys = {key: value for key, value in test_data.items() if key not in keys_to_ignore}

    # Print keys for debugging
    print("Keys in target_data_without_keys:", target_data_without_keys.keys())
    print("Keys in test_data_without_keys:", test_data_without_keys.keys())

    # Test each value
    # Conversion to string is required as float.nan != float.nan
    for key in target_data_without_keys.keys():
        if key in test_data_without_keys:
            assert str(test_data_without_keys[key]) == str(target_data_without_keys[key])
        else:
            print(f"Key {key} not found in test_data_without_keys")


if __name__ == "__main__":
    test_mixdesign_metadata_extraction()