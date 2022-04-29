
import unittest
import os
from pathlib import Path
from datetime import datetime

class test_generate_metadata_from_emodul_BAM(unittest.TestCase):
    
    # read specimen.dat file to make sure that headers work correctly
    def read_specimen_file(self, file):
        with open(file) as specimen:
            print ('test')
            #lines = specimen.readlines()
            # TODO: if checking line by line is needed it has to change
            for i, line in enumerate(specimen):
                if (specimen.name.__contains__('BA Los M VI-6')):
                    if ("Bediener" in line):
                        date_val = line.split('\t')
                        origin_date = date_val[9]
                        str_date = ''.join(origin_date.split())
                        str_date_time = datetime.strptime(str_date, '%d.%m.%Y%H:%M:%S')
                        self.assertIsInstance(str_date_time, datetime)
                        self.assertListEqual([str_date_time.year, str_date_time.month, str_date_time.day, 
                        str_date_time.hour, str_date_time.minute, str_date_time.second], 
                        [2014, 9, 11, 13, 22, 36])

                        zeit_val = date_val[7]
                        self.assertAlmostEqual(zeit_val, '17,580078')

                    if ("Probenbezeichnung " in line):
                        design_val = line.split(':')
                        str_design_val = ''.join(design_val[1].split())
                        self.assertEqual(str_design_val, 'BALosMVI-6')

                    if ("Masse" in line):
                        dimension = line.split(':')
                        dimension_val = ''.join(dimension[1].split())
                        self.assertEqual(dimension_val, '5367,6')
                    
                    if ("Durchmesser" in line):
                        diameter_val = line.split(':')
                        str_diameter_val = ''.join(diameter_val[1].split())
                        self.assertEqual(str_diameter_val, '98,7')

                    if ("Länge" in line):
                        length_val = line.split(':')
                        str_length_val = ''.join(length_val[1].split())
                        self.assertEqual(str_length_val, '300,3')

                    if ("Prüfer" in line):
                        examiner = line.split(':')
                        examiner_val = ''.join(examiner[1].split())
                        self.assertTrue(examiner_val.__eq__('Kh'))

                if (specimen.name.__contains__('Hüsken  E-Modul Probe 2-1')):
                    if ("Bediener" in line):
                        date_val = line.split('\t')
                        origin_date = date_val[9]
                        str_date = ''.join(origin_date.split())
                        str_date_time = datetime.strptime(str_date, '%d.%m.%Y%H:%M:%S')
                        self.assertIsInstance(str_date_time, datetime)
                        self.assertListEqual([str_date_time.year, str_date_time.month, str_date_time.day, 
                        str_date_time.hour, str_date_time.minute, str_date_time.second], 
                        [2014, 3, 18, 15, 46, 41])

                        zeit_val = date_val[7]
                        self.assertAlmostEqual(zeit_val, '56,78418')
                    if ("Probenbezeichnung " in line):
                        design_val = line.split(':')
                        str_design_val = ''.join(design_val[1].split())
                        self.assertEqual(str_design_val, 'HüskenE-ModulProbe2-1')

                    if ("Masse" in line):
                        dimension = line.split(':')
                        dimension_val = ''.join(dimension[1].split())
                        self.assertEqual(dimension_val, '5853')

                    if ("Durchmesser" in line):
                        diameter_val = line.split(':')
                        str_diameter_val = ''.join(diameter_val[1].split())
                        self.assertEqual(str_diameter_val, '100')

                    if ("Länge" in line):
                        length_val = line.split(':')
                        str_length_val = ''.join(length_val[1].split())
                        self.assertEqual(str_length_val, '300')

                    if ("Prüfer" in line):
                        examiner = line.split(':')
                        examiner_val = ''.join(examiner[1].split())
                        self.assertTrue(examiner_val.__eq__('Machura,Be'))

                
                # break in line 9 which is the end of the header
                if (i == 9):
                    break

    # read the directory content to find specimen.dat file
    def read_directory(self,path):
        filename = 'specimen.dat'
        with os.scandir(path) as items:
            for item in items:
                if item.name == filename:
                    self.read_specimen_file(item)


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
