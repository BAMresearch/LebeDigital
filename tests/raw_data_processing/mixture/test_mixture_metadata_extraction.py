import os
from lebedigital.raw_data_processing.mixture.mixdesign_metadata_extraction import extract_metadata_mixdesign
from pathlib import Path


def test_mixture_metadata_extraction():
    """Testing the mixture metadata extraction on a single example"""

    # setup paths and directories
    data_dir = 'test_data'
    data_path = Path(__file__).parent / data_dir
    input_file = '2014_08_05 Rezeptur_MI.xlsx'
    input_path = data_path / input_file
    output_file = data_path / 'processed_data.csv'

    # make sure that the path exists
    assert(os.path.exists(input_path))

    target_data = {"RawDataFile": "C:/develop/lebedigital-new/Lebedigital/usecases/MinimumWorkingExample/Demo_Workshop/Data/2014_08_05 Rezeptur_MI.xlsx",
    "MixingDate": "2014-06-30T12:00:00",
    "Lab": "BAM",
    "ID": "eac74cba-09fd-4155-b978-7eb570e6fcbe",
    "humanreadableID": "2014_08_05 Rezeptur_MI",
    "Cement1_Content": 330.0,
    "Cement1_Content_Unit": "kg/m^3",
    "Cement1_Density": 3.123,
    "Cement1_Density_Unit": "kg/dm^3",
    "Cement1_Type": "CEM I 42.5 R",
    "Water_Content": 175.0,
    "Water_Content_Unit": "kg/m^3",
    "Water_Density": 1.0,
    "Water_Density_Unit": "kg/dm^3",
    "WaterCementRatio": 0.5303030303030303,
    "Admixture1_Content": 5.61,
    "Admixture1_Content_Unit": "kg/m^3",
    "Admixture1_Density": 1.05,
    "Admixture1_Density_Unit": "kg/dm^3",
    "Admixture1_Type": "FM 595 BASF",
    "Aggregate1_Content": 1564.0,
    "Aggregate1_Content_Unit": "kg/m^3",
    "Aggregate1_Size": NaN,
    "Aggregate1_Size_Unit": NaN,
    "Aggregate1_Density": NaN,
    "Aggregate1_Density_Unit": "kg/dm^3",
    "Addition1_Content": 273.0,
    "Addition1_Content_Unit": "kg/m^3",
    "Addition1_Density": 2.74,
    "Addition1_Density_Unit": "kg/dm^3",
    "Addition1_Type": "Medenbach Kalksteinmehl"}

    # run extraction and getting a dictionary with metadata
    test_data = extract_metadata_mixdesign(input_path, None)

    # test each value
    # conversion to string is required as float.nan != float.nan
    for key in target_data.keys():
        assert(str(test_data[key]) == str(target_data[key]))
