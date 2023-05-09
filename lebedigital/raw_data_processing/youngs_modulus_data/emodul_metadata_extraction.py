import re
import json
import os
import sys
import pandas as pd
import yaml
from pathlib import Path
import argparse
import warnings
import uuid


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


def extract_metadata_emodulus(rawDataPath, specimen_file, mix_file):
    """Returns two dictionaries: one with extracted emodule-metadata and one with
    extracted specimen metadata.

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
    metadata_emodule : dict
        Return the dictionary with the extracted metadata for emodule
    metadata_specimen : dict
        Return the dictionary with specimen metadata
    
    """

    # create empty dictionary for metadata
    metadata_emodule = {}
    metadata_specimen = {}

    # read raw data file
    with open(str(rawDataPath)+'/'+str(specimen_file), encoding="utf8", errors='ignore') as data:


        # get metadata from file and location
        # get the name of the folder of the file
        folderName = os.path.basename(rawDataPath)

        # read data file line by line
        lines = data.readlines()

        # set software header - This data has no placeholder yet.
        metadata_emodule['software_specification'] = get_metadata_in_one_line(lines[0])[0] 

        # specific testing machine and software version this script is optimized for
        # - This data has no placeholder yet.
        assert metadata_emodule['software_specification'] == 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A' 

        # get empty lines (where start and end the header)
        emptyLineIndex = []
        for lineIndex in range(len(lines)):
            if len(lines[lineIndex]) == 1:
                emptyLineIndex.append(lineIndex)

        # service information of the experiment, it should be in between the first two empty lines
        serviceInformation = []
        for ind in range(emptyLineIndex[0]+1,emptyLineIndex[1]+4):
            serviceInformation.append(get_metadata_in_one_line(lines[ind]))



        ###########  D A T A   A B O U T    E X P E R I M E N T  #######

        # name of experiment is the folder name of the data file
        metadata_emodule['ExperimentName'] = folderName  # This data has no placeholder yet.

        # get experiment date and time
        date = serviceInformation[10][4]
        metadata_emodule['ExperimentDateTime'] = str(date)

        # operator name - This data has no placeholder yet.
        metadata_emodule['tester_name'] = serviceInformation[2][1] 

        # remarks - This data has no placeholder yet.
        metadata_emodule['remark'] = serviceInformation[4][1]

        # set experiment lab location to BAM
        metadata_emodule['Lab'] = 'BAM'

        # set Compression and Transducer Column
        metadata_emodule['CompressionColumn'] = 0
        metadata_emodule['AppliedLoad_Unit'] = "kN"
        metadata_emodule['TransducerColumn'] = [1,2,3]
        metadata_emodule['MeasuringGauge_Unit'] = "mm"  # Transducer messen eine Verschiebung.

        metadata_emodule['EModule_Value_Unit'] = None

        # set paths
        metadata_emodule['ProcessedFile'] = os.path.join('../usecases/MinimumWorkingExample/emodul/processed_data') # path to csv file with values extracted by emodul_generate_processed_data.py
        metadata_emodule['RawDataFile'] = os.path.join(rawDataPath,specimen_file) # path to specimen.dat
        try:
            with open(str(rawDataPath)+'/'+ str(mix_file), encoding="utf8", errors='ignore') as mix_data:
                lines = mix_data.readlines()
                lines = lines[0].strip()
                dataPath = Path(rawDataPath).parents[1]
                metadata_emodule['MixDataFile']= os.path.join(dataPath, "Mischungen", lines)
        except:
            metadata_emodule['MixDataFile'] = None


        ###########  D A T A   A B O U T    S P E C I M E N #######

        # name of specimen (human readable)
        metadata_emodule['SpecimenName'] = metadata_specimen['SpecimenName'] = serviceInformation[3][1] 

        # ID of specimen (machine readable)
        specimenID = str(uuid.uuid4())
        metadata_emodule['SpecimenID'] = metadata_specimen['SpecimenID'] = specimenID

        # set specimen age to 28 days
        metadata_emodule['SpecimenAge'] = 28.0
        metadata_emodule['SpecimenAge_Unit'] = 'day'

        # weight 
        metadata_specimen['SpecimenWeight'] = float(replace_comma(serviceInformation[5][1]))
        metadata_specimen['SpecimenWeight_Unit'] = 'g'

        # set size of specimen
        metadata_specimen['SpecimenDiameter'] = float(replace_comma(serviceInformation[6][1])) #diameter
        metadata_specimen['SpecimenDiameter_Unit'] = 'mm'
        metadata_specimen['SpecimenLength'] = float(replace_comma(serviceInformation[7][1])) #length
        metadata_specimen['SpecimenLength_Unit'] = 'mm'
        if metadata_specimen['SpecimenDiameter'] > metadata_specimen['SpecimenLength']:
            dir_name = metadata_emodule['ExperimentName']
            raise Exception(f'Diameter is larger then length, please fix the mistake in {dir_name}')

    return metadata_emodule, metadata_specimen


def emodul_metadata(rawDataPath, metaDataFile,specimenDataFile):
    """Creates two yaml files with extracted metadata, one for emodule and one
    for the specimen

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    metaDataFile : string
        Path to the output data file for emodule metadata
    specimenDataFile : string
        Path to the output data file for specimen metadata

    """

    # define file names for each data set
    mix_file = 'mix.dat'
    specimen_file = 'specimen.dat'

    # extracting the metadata
    metadata, specimen = extract_metadata_emodulus(rawDataPath, specimen_file, mix_file)
    
    # writing the metadata to yaml file
    with open(metaDataFile, 'w') as yamlFile:
        yaml.dump(metadata, yamlFile)
    with open(specimenDataFile, 'w') as yamlFile:
        yaml.dump(specimen, yamlFile)


def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from BAM e-module experiment')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata yaml
    parser.add_argument('-o', '--output', help='Path to extracted yaml files, one for emodule and one for specimen')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = '../../../usecases/MinimumWorkingExample/Data/E-modul/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4'
    if args.output == None:
        args.output = ['../../../usecases/MinimumWorkingExample/emodul/metadata_yaml_files/testMetaData.yaml','../../../usecases/MinimumWorkingExample/emodul/metadata_yaml_files/testSpecimenData.yaml']

    # run extraction and write metadata file
    emodul_metadata(args.input, args.output[0],args.output[1])


if __name__ == "__main__":
    main()


