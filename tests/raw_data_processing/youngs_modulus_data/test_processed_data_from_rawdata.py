import unittest
import os

from lebedigital.raw_data_processing.youngs_modulus_data.emodul_generate_processed_data import \
    processed_data_from_rawdata

class TestProcessedDataFromRawData(unittest.TestCase):

    def test_emodul_generate_processed_data(self):
        """
        Testing the processed data generation on a single experiment
        """

        rawdata_location = '../usecases/MinimumWorkingExample/Data/E-modul/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4'
        processeddat_location = './BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4' + '.csv'

        processed_data_from_rawdata(locationOfRawData=rawdata_location, locationOfProcessedData=processeddat_location)

        self.assertTrue(os.path.exists('./BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.csv'))

