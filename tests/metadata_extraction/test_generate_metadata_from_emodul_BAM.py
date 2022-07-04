
from email.contentmanager import raw_data_manager
import unittest
import os
from pathlib import Path
from datetime import datetime

from lebedigital.raw_data_processing.metadata_extraction.emodul_metadata_extraction import eModul_metadata

class test_generate_metadata_from_emodul_BAM(unittest.TestCase):

    # read the directory content to find specimen.dat file
    def read_directory(self,path):
        filename = 'specimen.dat'
        metadata = []
        # call emodul_metadata function from the script
        metadata.append(eModul_metadata(path, filename, ''))
        return metadata

    def test_BA_Los_M_VI_6 (self):
        print('test ba los m vi 6')
        baseDir = Path(__file__).resolve().parents[1]
        dataPath = os.path.join(baseDir, 'usecases', 'MinimumWorkingExample','Data', 'E-modul')

        # walk through the directories
        with os.scandir(dataPath) as ls:
            sub_dir = 'BA Los M VI-6'
            inputParam = os.path.join(dataPath, sub_dir)
            metadata = self.read_directory(inputParam)

            # read specimen.dat file to make sure that headers work correctly
            date_val = metadata[0]['operatorTimestamp']
            str_date_time = datetime.strptime(date_val, '%d.%m.%Y %H:%M:%S')
            self.assertIsInstance(str_date_time, datetime)
            self.assertListEqual([str_date_time.year, str_date_time.month, str_date_time.day, 
            str_date_time.hour, str_date_time.minute, str_date_time.second], 
            [2014, 9, 11, 13, 22, 36])

            returned_metadata = metadata[0]

            zeit_val = returned_metadata['operatorTime']
            self.assertAlmostEqual(zeit_val, '17,580078')

            design_val = returned_metadata['experimentName']
            self.assertEqual(design_val, 'BA Los M VI-6')

            dimension = returned_metadata['weight']
            self.assertEqual(dimension, '5367,6')
            
            diameter_val = returned_metadata['diameter']
            self.assertEqual(diameter_val, '98,7')

            length_val = returned_metadata['length']
            self.assertEqual(length_val, '300,3')

            examiner = returned_metadata['tester']
            self.assertTrue(examiner.__eq__('Kh'))

    def test_Hüsken_EModul_Probe (self):
        print('test Hüsken  E-Modul Probe 2-1')
        baseDir = Path(__file__).resolve().parents[1]
        dataPath = os.path.join(baseDir, 'usecases', 'MinimumWorkingExample','Data', 'E-modul')
        file_paths = []

        # walk through the directories
        with os.scandir(dataPath) as ls:
            sub_dir = 'Hüsken  E-Modul Probe 2-1'
            inputParam = os.path.join(dataPath, sub_dir)
            metadata = self.read_directory(inputParam)

            # read specimen.dat file to make sure that headers work correctly
            date_val = metadata[0]['operatorTimestamp']
            str_date_time = datetime.strptime(date_val, '%d.%m.%Y %H:%M:%S')
            self.assertIsInstance(str_date_time, datetime)
            self.assertListEqual([str_date_time.year, str_date_time.month, str_date_time.day, 
            str_date_time.hour, str_date_time.minute, str_date_time.second], 
            [2014, 3, 18, 15, 46, 41])
            returned_metadata = metadata[0]

            zeit_val = returned_metadata['operatorTime']
            self.assertAlmostEqual(zeit_val, '56,78418')

            design_val = returned_metadata['experimentName']
            self.assertEqual(design_val, 'Hüsken  E-Modul Probe 2-1')

            dimension = returned_metadata['weight']
            self.assertEqual(dimension, '5853')

            diameter_val = returned_metadata['diameter']
            self.assertEqual(diameter_val, '100')
            
            length_val = returned_metadata['length']
            self.assertEqual(length_val, '300')

            examiner = returned_metadata['tester']
            self.assertTrue(examiner.__eq__('Machura, Be'))

if __name__ == '__main__':
     unittest.main()
