import pytest
import os

from lebedigital.raw_data_processing.processed_data_generation.emodul_generate_processed_data import \
    processed_data_from_rawdata


def test_emodul_generate_processed_data():
    """
    Testing the processed data generation on a single experiment
    """

    rawdata_location = '../usecases/MinimumWorkingExample/Data/E-modul/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4'
    processeddat_location = './BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4' + '.csv'

    processed_data_from_rawdata(locationOfRawData=rawdata_location, locationOfProcessedData=processeddat_location)

    assert os.path.exists('./BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.csv') == True