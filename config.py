import os


class Config:
    """
        This class defines some general configs

        Attributes
        ----------
        TEST_MIXTURE_METADATA_EXTRACTION_EXPECTED_VALUE: str
            get the expected value to compare with the actual one
        """

    TEST_MIXTURE_METADATA_EXTRACTION_EXPECTED_VALUE = "{'operator_date': '2014-06-30', 'tester_name': 'Werner', 'specimen_name': 'BA-Losert MI', 'cement--QuantityInMix': 330.0, 'cement--BulkDensity': 3.123, 'cement--Volume': 105.7, 'cement--Annotation': 'CEM I 42.5 R', 'water_total--QuantityInMix': 175.0, 'water_total--BulkDensity': 1.0, 'water_total--Volume': 175.0, 'water_cement_ratio': 0.5303030303030303, 'water_effective--QuantityInMix': nan, 'water_effective--BulkDensity': nan, 'water_effective--Volume': nan, 'air_content--QuantityInMix': 0.0, 'air_content--BulkDensity': 0.0, 'air_content--Volume': 20.0, 'addition1--QuantityInMix': 273.0, 'addition1--BulkDensity': 2.74, 'addition1--Volume': 99.6, 'addition1--Annotation': 'Medenbach - Kalksteinmehl', 'admixture--QuantityInMix': 5.61, 'admixture--BulkDensity': 1.05, 'admixture--Volume': 5.3, 'admixture--Annotation': 'FM 595 BASF', 'aggregate--QuantityInMix': 1564.0, 'aggregate--BulkDensity': nan, 'aggregate--Volume': 594.4000000000001}"

    TEST_METADATA_EXTRACTION_EXPECTED_VALUE = {'experimentName': 'BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4',
                                               'software_specification': 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A',
                                               'operator_timestamp': '13:25:39',
                                               'operator_date': '01.09.2014',
                                               'tester_name': 'Kh',
                                               'specimen_name': 'BA-Losert E-Modul 28d v. 04.08.14 Probe 4',
                                               'remark': 'Kraftgeregelt 3,9 kN/s',
                                               'weight': 5342.0,
                                               'weight_unit': 'g',
                                               'diameter': 98.6,
                                               'length': 300.3,
                                               'length_unit': 'mm',
                                               'mix_file': '2014_08_05 Rezeptur_MI.xlsx'}
