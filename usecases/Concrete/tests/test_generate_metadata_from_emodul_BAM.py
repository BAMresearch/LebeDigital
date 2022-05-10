
import unittest
import os
from pathlib import Path
from datetime import datetime
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from knowledgeGraph.emodul.emodul_metadata_extraction import eModul_metadata as emodul

class test_generate_metadata_from_emodul_BAM(unittest.TestCase):
    
    # read specimen.dat file to make sure that headers work correctly
    def read_specimen_file(self, metadata):
        # separate between two headers of selected files
        if (metadata[0][0]['experimentName'] == 'BA Los M VI-6'):
            date_val = metadata[0][2]['Bediener Information']['Zeitpunkt']
            str_date_time = datetime.strptime(date_val, '%d.%m.%Y %H:%M:%S')
            self.assertIsInstance(str_date_time, datetime)
            self.assertListEqual([str_date_time.year, str_date_time.month, str_date_time.day, 
            str_date_time.hour, str_date_time.minute, str_date_time.second], 
            [2014, 9, 11, 13, 22, 36])

            zeit_val = metadata[0][2]['Bediener Information']['Zeit']
            self.assertAlmostEqual(zeit_val, '17,580078')

            design_val = metadata[0][2]['Bediener Information']['Probenbezeichnung']
            self.assertEqual(design_val, 'BA Los M VI-6')

            dimension = metadata[0][2]['Bediener Information']['Masse']
            self.assertEqual(dimension, '5367,6')
            
            diameter_val = metadata[0][2]['Bediener Information']['Durchmesser']
            self.assertEqual(diameter_val, '98,7')

            length_val = metadata[0][2]['Bediener Information']['Lnge']
            self.assertEqual(length_val, '300,3')

            examiner = metadata[0][2]['Bediener Information']['Prfer']
            self.assertTrue(examiner.__eq__('Kh'))

        if (metadata[0][0]['experimentName'] == 'Hüsken  E-Modul Probe 2-1'):
            date_val = metadata[0][2]['Bediener Information']['Zeitpunkt']
            str_date_time = datetime.strptime(date_val, '%d.%m.%Y %H:%M:%S')
            self.assertIsInstance(str_date_time, datetime)
            self.assertListEqual([str_date_time.year, str_date_time.month, str_date_time.day, 
            str_date_time.hour, str_date_time.minute, str_date_time.second], 
            [2014, 3, 18, 15, 46, 41])

            zeit_val = metadata[0][2]['Bediener Information']['Zeit']
            self.assertAlmostEqual(zeit_val, '56,78418')

            design_val = metadata[0][2]['Bediener Information']['Probenbezeichnung']
            self.assertEqual(design_val, 'Hsken  E-Modul Probe 2-1')

            dimension = metadata[0][2]['Bediener Information']['Masse']
            self.assertEqual(dimension, '5853')

            diameter_val = metadata[0][2]['Bediener Information']['Durchmesser']
            self.assertEqual(diameter_val, '100')
            
            length_val = metadata[0][2]['Bediener Information']['Lnge']
            self.assertEqual(length_val, '300')

            examiner = metadata[0][2]['Bediener Information']['Prfer']
            self.assertTrue(examiner.__eq__('Machura, Be'))


    # read the directory content to find specimen.dat file
    def read_directory(self,path):
        filename = 'specimen.dat'
        metadata = []
        # call emodul_metadata function from the script
        metadata.append(emodul(path, filename))
        self.read_specimen_file(metadata)


    def test_generate_metadata_yaml_file_from_emodul_BAM(self):
        baseDir = Path(__file__).resolve().parents[1]
        dataPath = os.path.join(baseDir, 'Data', 'E-modul')
        file_paths = []

        # walk through the directories
        with os.scandir(dataPath) as ls:
            sub_dirs = ['BA Los M VI-6', 'Hüsken  E-Modul Probe 2-1']
            for item in sub_dirs:
                inputParam = os.path.join(dataPath, item)
                self.read_directory(inputParam)
        return file_paths

if __name__ == '__main__':
    unittest.main()
