import re
import json
import os
import sys
import pandas as pd
import yaml
from pathlib import Path
import argparse
import warnings


# the function read each line and return metadata as key and value
def get_metadata_in_one_line(line):
    s = re.sub('\t+', '\t', line)
    s = s.replace('\n','\t')
    result = s.split('\t')[:-1]
    return result


# function to convert german formatting to english
def replace_comma(string):
    string = string.replace(',', '.')
    return string


def extract_metadata_emodulus(rawDataPath,specimen_file,mix_file):
    """Return a dictionary with extracted metadata

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    specimen_file : string
        Name of the file containing most meta data
    mix_file : string
        File name containing the name of the file with the mix data

    Returns
    -------
    metadata : dict
        Return the dictionary with the extracted metadata
    """

    # create empty dictionary for metadata
    metadata = {}

    # read raw data file
    with open(str(rawDataPath)+'/'+str(specimen_file), encoding="utf8", errors='ignore') as data:


        # get metadata from file and location
        # get the name of the folder of the file
        folderName = os.path.basename(rawDataPath)

        # read data file line by line
        lines = data.readlines()

        # set software header - This data has no placeholder yet.
        metadata['software_specification'] = get_metadata_in_one_line(lines[0])[0] 

        # specific testing machine and software version this script is optimized for
        # - This data has no placeholder yet.
        assert metadata['software_specification'] == 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A' 

        # get empty lines (where start and end the header)
        emptyLineIndex = []
        for lineIndex in range(len(lines)):
            if len(lines[lineIndex]) == 1:
                emptyLineIndex.append(lineIndex)

        # service information of the experiment, it should be in between the first two empty lines
        serviceInformation = []
        for ind in range(emptyLineIndex[0]+1,emptyLineIndex[1]+2):
            serviceInformation.append(get_metadata_in_one_line(lines[ind]))



        ###########  D A T A   A B O U T    E X P E R I M E N T  #######

        # name of experiment is the folder name of the data file
        metadata['experimentName'] = folderName  # This data has no placeholder yet.

        # get experiment date and time
        date, time = serviceInformation[0][4].split(' ')
        metadata['ExperimentTime'] = str(time) # operator_timestamp
        metadata['ExperimentDate'] = str(date) # operator_date

        # get measurement duration  
        metadata['MeasurementDuration'] = float(replace_comma(serviceInformation[10][2]))

        # operator name - This data has no placeholder yet.
        metadata['tester_name'] = serviceInformation[2][1] 

        # remarks - This data has no placeholder yet.
        metadata['remark'] = serviceInformation[4][1]

        # set experiment lab location to BAM
        metadata['Lab'] = 'BAM'

        # set Compression and Transducer Column
        metadata['CompressionColumn'] = 0
        metadata['TransducerColumn'] = [1,2,3]

        # set paths
        metadata['ProcessedFile'] = os.path.join('../usecases/MinimumWorkingExample/emodul/processed_data') # path to csv file with values extracted by emodul_generate_processed_data.py
        metadata['RawDataFile'] = os.path.join(rawDataPath,specimen_file) # path to specimen.dat
        try:
            with open(str(rawDataPath)+'/'+ str(mix_file), encoding="utf8", errors='ignore') as mix_data:
                lines = mix_data.readlines()
                lines = lines[0].strip()
                dataPath = Path(rawDataPath).parents[1]
                metadata['MixDataFile']= os.path.join(dataPath, "Mischungen", lines)
        except:
            metadata['MixDataFile'] = None


        ###########  D A T A   A B O U T    S P E C I M E N #######

        # name of specimen - This data has no placeholder yet.
        metadata['specimen_name'] = serviceInformation[3][1] 

        # set specimen age to 28 days
        metadata['SpecimenAge'] = 28.0

        # weight - This data has no placeholder yet.
        weight = float(replace_comma(serviceInformation[5][1]))
        metadata['weight'] = weight

        # set size of specimen
        metadata['SpecimenDiameter'] = float(replace_comma(serviceInformation[6][1])) #diameter
        metadata['SpecimenLength'] = float(replace_comma(serviceInformation[7][1])) #length
        if metadata['SpecimenDiameter'] > metadata['SpecimenLength']:
            dir_name = metadata['experimentName']
            raise Exception(f'Diameter is larger then length, please fix the mistake in {dir_name}')

        # set units, assuming length and diameter have the same - This data has no placeholder yet.
        if metadata['SpecimenLength'] > 100:
            metadata['length_unit'] = 'mm'
        elif metadata['SpecimenLength'] < 1:
            metadata['length_unit'] = 'm'
        else:
            raise Exception('Unexpected value of length')

    return metadata


def emodul_metadata(rawDataPath, metaDataFile):
    """Creates a yaml file with extracted metadata

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    metaDataFile : string
        Path to the output data file
    """

    # define file names for each data set
    mix_file = 'mix.dat'
    specimen_file = 'specimen.dat'

    # extracting the metadata
    metadata = extract_metadata_emodulus(rawDataPath, specimen_file, mix_file)

    # writing the metadata to yaml file
    with open(metaDataFile, 'w') as yamlFile:
        yaml.dump(metadata, yamlFile)


def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from BAM e-module experiment')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata yaml
    parser.add_argument('-o', '--output', help='Path to extracted meta data yaml file')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = '../../../usecases/MinimumWorkingExample/Data/E-modul/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4'
    if args.output == None:
        args.output = '../../../usecases/MinimumWorkingExample/Data/testMetaData.yaml'
        #args.output = '../../../usecases/MinimumWorkingExample/emodul/metadata_yaml_files/testMetaData.yaml'

    # run extraction and write metadata file
    emodul_metadata(args.input, args.output)


if __name__ == "__main__":
    main()


