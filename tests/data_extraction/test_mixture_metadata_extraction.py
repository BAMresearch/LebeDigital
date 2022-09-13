import pytest
from lebedigital.raw_data_processing.metadata_extraction \
    .mixture_metadata_extraction import extraction


def test_extraction():
    """Tesing the mixture metadata extraction on a single example"""

    # setting up the test example
    input = '../usecases/MinimumWorkingExample/Data/Mischungen/2014_08_05 Rezeptur_MIII.xlsx'
    target_data = {'Cement -- Quantity in mix': 330, 'Cement -- Specific gravity': 3.123,
                   'Cement -- Annotation': 'CEM I 42,5 R', 'Water (total) -- Quantity in mix':
                    175, 'Water (total) -- Specific gravity': 1, 'Water (total) -- Annotation': 
                    '---', 'Water (effective) -- Quantity in mix': 
                    '---', 'Water (effective) -- Specific gravity': 
                    '---', 'Water (effective) -- Annotation': '---', 'Air -- Quantity in mix': 
                    '---', 'Air -- Specific gravity': '---', 'Air -- Annotation': '---', 
                    'Supplementary cementious materials Kalksteinmehl -- Quantity in mix': 147, 
                    'Supplementary cementious materials Kalksteinmehl -- Specific gravity': 2.74,
                    'Supplementary cementious materials Kalksteinmehl -- Annotation': 
                    'Medenbach', 'Admixture -- Quantity in mix': 4.95, 
                    'Admixture -- Specific gravity': 1.05, 'Admixture -- Annotation': 
                    'FM 595 BASF', 'Total -- Quantity in mix': 656.95, 
                    'Total -- Specific gravity': '---', 'Total -- Annotation': '---', 
                    'Zuschlag (total) -- Quantity in mix': 1685, 
                    'Zuschlag (total) -- Specific gravity': '---', 
                    'Zuschlag (total) -- Annotation': '---', '0 / 0,3 -- Quantity in mix': 0, 
                    '0 / 0,3 -- Specific gravity': 2.65, '0 / 0,3 -- Annotation': 'Quarz', 
                    '0,1 / 0,5 -- Quantity in mix': 337, '0,1 / 0,5 -- Specific gravity': 2.63, 
                    '0,1 / 0,5 -- Annotation': 'Okrilla  ', '0,5 / 1,0 -- Quantity in mix': 194, 
                    '0,5 / 1,0 -- Specific gravity': 2.63,
                    '0,5 / 1,0 -- Annotation': 'Okrilla  ', '1,0 / 2,0 -- Quantity in mix': 194, 
                    '1,0 / 2,0 -- Specific gravity': 2.63, '1,0 / 2,0 -- Annotation': 'Okrilla  ', 
                    '2,0 / 4,0 -- Quantity in mix': 126, '2,0 / 4,0 -- Specific gravity': 2.63, 
                    '2,0 / 4,0 -- Annotation': 'Okrilla  ', '4,0 / 8,0 -- Quantity in mix': 126, 
                    '4,0 / 8,0 -- Specific gravity': 2.63, '4,0 / 8,0 -- Annotation': 'Okrilla  ', 
                    '8,0 / 16,0 -- Quantity in mix': 708, '8,0 / 16,0 -- Specific gravity': 2.63, 
                    '8,0 / 16,0 -- Annotation': 'Okrilla  ', 'nan -- Quantity in mix': '---', 
                    'nan -- Specific gravity': '---', 'nan -- Annotation': '---', 
                    'Fresh concrete -- Quantity in mix': 2341.95, 
                    'Fresh concrete -- Specific gravity': '---', 'Fresh concrete -- Annotation': 
                    '---', 'Mehlkornanteil -- Quantity in mix': '---', 
                    'Mehlkornanteil -- Specific gravity': '---', 'Mehlkornanteil -- Annotation': 
                    '---', 'Mörtelanteil -- Quantity in mix': '---', 
                    'Mörtelanteil -- Specific gravity': '---', 'Mörtelanteil -- Annotation': '---', 
                    'Date': '2014-06-30', 'Editor': 'Werner', 'Requester': 'BAM 7.0', 
                    'Project number': 'BA-Losert', 'Specimen': 'BA-Losert MIII', 
                    'Betonsorte u -festigkeitsklasse': '---', 'Wasserzementwert': '0,53 w/b = 0,37', 
                    'Konsistenz': 'SVB', 'Sieblinie n. DIN 1045': 'SVB', 'Körnungsziffer': 3.88, 
                    'Vol Leim/Zuschlag': 56.00624024960999}
    # run extraction and getting a dictionary with metadata
    test_data = extraction(input, None)

    # checking if result is correct
    assert test_data == target_data
