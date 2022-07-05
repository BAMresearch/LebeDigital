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

        # name of experiment is the folder name of the data file
        metadata['experimentName'] = folderName


        # read data file line by line
        lines = data.readlines()

        # set software header
        metadata['software_specification'] = get_metadata_in_one_line(lines[0])[0]

        # specific testing machine and software version this script is optimized for
        assert metadata['software_specification'] == 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A'

        # get empty lines (where start and end the header)
        emptyLineIndex = []
        for lineIndex in range(len(lines)):
            if len(lines[lineIndex]) == 1:
                emptyLineIndex.append(lineIndex)

        # service information of the experiment, it should be in between the first two empty lines
        serviceInformation = []
        for ind in range(emptyLineIndex[0]+1,emptyLineIndex[1]):
            serviceInformation.append(get_metadata_in_one_line(lines[ind]))

        # get date and time
        date, time = serviceInformation[0][4].split(' ')

        metadata['operator_timestamp'] = str(time)
        metadata['operator_date'] = str(date)

        # operator name
        metadata['tester_name'] = serviceInformation[2][1]

        # name of specimen
        metadata['specimen_name'] = serviceInformation[3][1]

        # remarks
        metadata['remark'] = serviceInformation[4][1]

        # weight
        weight = float(replace_comma(serviceInformation[5][1]))
        metadata['weight'] = weight

        # set weight unit, only a rough approximation
        # should be included in the input data in the long run
        if metadata['weight'] > 1000:
            metadata['weight_unit'] = 'g'
        elif metadata['weight'] < 10:
            metadata['weight_unit'] = 'kg'
        else:
            raise Exception('Unexpected value of weight')

        # set size of specimen
        metadata['diameter'] = float(replace_comma(serviceInformation[6][1]))
        metadata['length'] = float(replace_comma(serviceInformation[7][1]))

        if metadata['diameter'] > metadata['length']:
            dir_name = metadata['experimentName']
            raise Exception(f'Diameter is larger then length, please fix the mistake in {dir_name}')

        # set units, assuming length and diameter have the same
        if metadata['length'] > 100:
            metadata['length_unit'] = 'mm'
        elif metadata['length'] < 1:
            metadata['length_unit'] = 'm'
        else:
            raise Exception('Unexpected value of length')

        try:
            with open(str(rawDataPath)+'/'+str(mix_file), encoding="utf8", errors='ignore') as mix_data:
                lines = mix_data.readlines()
                metadata['mix_file'] = lines[0].strip()
        except:
            metadata['mix_file'] = None

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
        args.output = '../../../usecases/MinimumWorkingExample/emodul/metadata_yaml_files/testMetaData.yaml'

    # run extraction and write metadata file
    emodul_metadata(args.input, args.output)


if __name__ == "__main__":
    main()


