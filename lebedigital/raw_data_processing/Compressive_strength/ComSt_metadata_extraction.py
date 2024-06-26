import re
import json
import os
import argparse
import warnings
import uuid
import datetime
import pandas as pd


# the function read each line and return metadata as key and value
def get_metadata_in_one_line(line):
    s = re.sub('\t+', '\t', line)
    s = s.replace('\n', '\t')
    result = s.split('\t')[:-1]
    return result


# function to convert german formatting to english
def replace_comma(string):
    string = string.replace(',', '.')
    return string


def extract_metadata_ComSt(rawDataPath, specimen_file='specimen.dat', mix_file='mix.dat'):
    """Returns two dictionaries: one with extracted emodule-metadata and one with
    extracted specimen metadata.

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    specimen_file : string
        Name of the file containing most metadata
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
    metadata_ComSt = {}
    metadata_specimen_ComSt = {}

    # read raw data file
    #with open(str(rawDataPath), encoding="utf8", errors='ignore') as data:
    with open(str(rawDataPath) + '/' + str(specimen_file), encoding="utf8", errors='ignore') as data:

        # get metadata from file and location
        # get the name of the folder of the file
        folderName = os.path.basename(rawDataPath)

        # read data file line by line
        lines = data.readlines()

        # set software header - This data has no placeholder yet.
        software_header = get_metadata_in_one_line(lines[0])[0]

        # specific testing machine and software version this script is optimized for
        # - This data has no placeholder yet.
        # assert software_header == 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A'

        # get empty lines (where start and end the header)
        emptyLineIndex = []
        for lineIndex in range(len(lines)):
            if len(lines[lineIndex]) == 1:
                emptyLineIndex.append(lineIndex)

        # service information of the experiment, it should be in between the first two empty lines
        serviceInformation = []
        for ind in range(emptyLineIndex[0] + 1, emptyLineIndex[1] + 4):
            serviceInformation.append(get_metadata_in_one_line(lines[ind]))

        ###########  D A T A   A B O U T    E X P E R I M E N T  #######

        # humanreadable ID = name of experiment is the folder name of the data file
        metadata_ComSt['humanreadableID'] = folderName

        # ID of this experiment
        ComStID = str(uuid.uuid4())
        metadata_ComSt['ID'] = ComStID

        # get experiment date and time in protege format YYYY-MM-DDTHH:mm:SS
        date = serviceInformation[11][4]  # datetime.datetime.strptime(,'%d.%m.%y')
        date_only = datetime.datetime.strptime(date.split(" ")[0], '%d.%m.%Y')
        date_protegeformat = date_only.strftime('%Y-%m-%d') + "T" + date.split(" ")[1]
        metadata_ComSt['ExperimentDate'] = str(date_protegeformat)

        # operator name - This data has no placeholder yet.
        # metadata_ComSt['tester_name'] = serviceInformation[2][1]

        # remarks - This data has no placeholder yet.
        # metadata_ComSt['remark'] = serviceInformation[4][1]

        # set experiment lab location to BAM
        metadata_ComSt['Lab'] = 'BAM'

        # set Compression and Transducer Column
        metadata_ComSt['CompressionColumn'] = [4]
        metadata_ComSt['CompressionForce_Unit'] = "kN"
        metadata_ComSt['TransducerColumn'] = [5]
        metadata_ComSt['Extensometer_Unit'] = "mm"  # Transducer messen eine Verschiebung.


        # name of specimen (humanreadable)
        metadata_specimen_ComSt['humanreadableID'] = folderName
        # set size of specimen
        metadata_specimen_ComSt['SpecimenDiameter'] = float(replace_comma(serviceInformation[5][1]))  # diameter
        metadata_specimen_ComSt['SpecimenDiameter_Unit'] = 'mm'
        metadata_specimen_ComSt['SpecimenHeight'] = float(replace_comma(serviceInformation[6][1]))  # Height
        metadata_specimen_ComSt['SpecimenHeight_Unit'] = 'mm'

        # set the length
        metadata_specimen_ComSt['SpecimenLength'] = float(replace_comma(serviceInformation[7][1]))  # Length
        metadata_specimen_ComSt['SpecimenLength_Unit'] = 'mm'

        # weight
        metadata_specimen_ComSt['SpecimenMass'] = float(replace_comma(serviceInformation[8][1]))
        metadata_specimen_ComSt['SpecimenMass_Unit'] = 'g'

        # path to specimen.dat
        try:
            #with open(str(rawDataPath), encoding="utf8", errors='ignore') as mix_data:
            with open(str(rawDataPath) + '/' + str(mix_file), encoding="utf8", errors='ignore') as mix_data:
            #with open(path_to_json) as mix_data:
                lines = mix_data.readlines()
                lines = lines[0].strip()

        except:
            # metadata_emodule['MixDataFile'] = None
            raise Exception("No mixdesign json-file found!")

        # ID of this specimen
        #specimenID = str(uuid.uuid4())
        metadata_ComSt['specimenID'] = metadata_specimen_ComSt['ID'] = ComStID
        # save Mixdesign ID to specimen metadata
        try:
            with open("../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/" + os.path.splitext(lines)[0] + ".json", "r", encoding="utf8", errors='ignore') as mixjson:  # change location of where
            #with open(path_to_json, "r") as mixjson:
                mixdesign = json.load(mixjson)
                mixtureID = mixdesign['ID']
                metadata_specimen_ComSt['MixtureID'] = mixtureID
            # Extract mixing date
            mixing_date_str = mixdesign['MixingDate']
            mixing_date = datetime.datetime.strptime(mixing_date_str, '%Y-%m-%dT%H:%M:%S')
        except AttributeError:
            raise Exception("No mixdesign json-file found! Can't import the ID and save it to the output!")

        # Calculate specimen age
        if 'ExperimentDate' in metadata_ComSt and mixing_date:
            # Convert dates to midnight
            experiment_date = datetime.datetime.strptime(metadata_ComSt['ExperimentDate'], '%Y-%m-%dT%H:%M:%S').replace(
                hour=0, minute=0, second=0)
            mixing_date = mixing_date.replace(hour=0, minute=0, second=0)
            # Calculate age
            specimen_age = (experiment_date - mixing_date).days
            metadata_ComSt['SpecimenAge'] = specimen_age
            metadata_ComSt['SpecimenAge_Unit'] = 'day'

        # set shape
        if 'SpecimenHeight' in metadata_specimen_ComSt:
            metadata_specimen_ComSt['SpecimenShape'] = 'Cube'

        # set paths
        metadata_ComSt['ProcessedFile'] = os.path.join('../../../usecases/MinimumWorkingExample/Druckfestigkeit/processeddata')  # path to csv file with values extracted by ComSt_generate_processed_data.py
        metadata_ComSt['RawDataFile'] = os.path.join(rawDataPath, specimen_file).replace('\\', '/')

        try:
            my_data = pd.read_csv(metadata_ComSt['ProcessedFile'])
            # print(my_data)
            min_force = my_data['Force [kN]'].min()
            #diameter = 100.0
            #height = 100.3
            area = metadata_specimen_ComSt['SpecimenDiameter'] * metadata_specimen_ComSt['SpecimenDiameter']
            normalValue = (min_force * -1)
            CompressiveStrength = (normalValue / area) * 1000
            metadata_ComSt['CompressiveStrength'] = CompressiveStrength
        except:

            raise Exception("No processed_file found!")

        metadata_ComSt['CompressiveStrength_Unit'] = "GPa"

    return metadata_ComSt, metadata_specimen_ComSt


def ComSt_metadata(rawDataPath, metaDataFile, specimenDataFile):
    """Creates two json files with extracted metadata, one for ComSt and one
    for the specimen

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    metaDataFile : string
        Path to the output data file for ComSt metadata
    specimenDataFile : string
        Path to the output data file for specimen metadata

    """

    # define file name for data set
    mix_file = 'mix.dat'
    specimen_file = 'specimen.dat'

    # extracting the metadata
    metadata, specimen = extract_metadata_ComSt(rawDataPath, specimen_file, mix_file)

    with open(metaDataFile, 'w') as jsonFile:
        json.dump(metadata, jsonFile, sort_keys=False, ensure_ascii=False, indent=4)
    with open(specimenDataFile, 'w') as jsonFile:
        json.dump(specimen, jsonFile, sort_keys=False, ensure_ascii=False, indent=4)


def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from BAM ComSt experiment')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata json
    parser.add_argument('-o', '--output', help='Path to extracted json files, one for ComSt and one for specimen')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = '../../../usecases/MinimumWorkingExample/Data/Druckfestigkeit_BAM/20240220_7188_M01/Druckfestigkeiten_7Tage/20240220_7188_M01_W01'
    if args.output == None:
        args.output = ['../../../usecases/MinimumWorkingExample/Druckfestigkeit/metadata_json_files/20240220_7188_M01_W01_ComSt.json',
                       '../../../usecases/MinimumWorkingExample/Druckfestigkeit/metadata_json_files/20240220_7188_M01_W01_ComStSpecimen.json']

    # run extraction and write metadata file
    ComSt_metadata(args.input, args.output[0], args.output[1])


if __name__ == "__main__":
    main()

