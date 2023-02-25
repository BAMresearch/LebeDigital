import os
from lebedigital.raw_data_processing.mixture.mixture_metadata_extraction import extract_metadata_mixture


def test_mixture_metadata_extraction():
    """Testing the mixture metadata extraction on a single example"""

    # setting up the test example input
    input_path = '../usecases/MinimumWorkingExample/Data/Mischungen/2014_08_05 Rezeptur_MI.xlsx'

    # make sure that the path exists
    assert(os.path.exists(input_path))

    target_data = {'operator_date': '2014-06-30',
                   'tester_name': 'Werner',
                   'specimen_name': 'BA-Losert MI',
                   'cement--QuantityInMix': 330.0,
                   'cement--BulkDensity': 3.123,
                   'cement--Volume': 105.7,
                   'cement--Annotation': 'CEM I 42.5 R',
                   'water_total--QuantityInMix': 175.0,
                   'water_total--BulkDensity': 1.0,
                   'water_total--Volume': 175.0,
                   'water_cement_ratio': 0.5303030303030303,
                   'water_effective--QuantityInMix': float('nan'),
                   'water_effective--BulkDensity': float('nan'),
                   'water_effective--Volume': float('nan'),
                   'air_content--QuantityInMix': 0.0,
                   'air_content--BulkDensity': 0.0,
                   'air_content--Volume': 20.0,
                   'addition1--QuantityInMix': 273.0,
                   'addition1--BulkDensity': 2.74,
                   'addition1--Volume': 99.6,
                   'addition1--Annotation': 'Medenbach - Kalksteinmehl',
                   'admixture--QuantityInMix': 5.61,
                   'admixture--BulkDensity': 1.05,
                   'admixture--Volume': 5.3,
                   'admixture--Annotation': 'FM 595 BASF',
                   'aggregate--QuantityInMix': 1564.0,
                   'aggregate--BulkDensity': float('nan'),
                   'aggregate--Volume': 594.4000000000001}

    # run extraction and getting a dictionary with metadata
    test_data = extract_metadata_mixture(input_path, None)

    # test each value
    # conversion to string is required as float.nan != float.nan
    for key in target_data.keys():
        assert(str(test_data[key]) == str(target_data[key]))
