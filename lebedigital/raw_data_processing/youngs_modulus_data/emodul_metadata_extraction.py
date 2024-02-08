import re
import json
import os
import argparse
import warnings
import uuid
import datetime


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


def extract_metadata_emodulus(rawDataPath, specimen_file='specimen.dat', mix_file='mix.dat'):

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
    metadata_emodule = {}
    metadata_specimen = {}

    # read raw data file
    #with open(str(rawDataPath), encoding="utf8", errors='ignore') as data:
    with open(str(rawDataPath)+'/'+str(specimen_file), encoding="utf8", errors='ignore') as data:


        # get metadata from file and location
        # get the name of the folder of the file
        folderName = os.path.basename(rawDataPath)

        # read data file line by line
        lines = data.readlines()

        # set software header - This data has no placeholder yet.
        software_header = get_metadata_in_one_line(lines[0])[0]

        # specific testing machine and software version this script is optimized for
        # - This data has no placeholder yet.
        #assert software_header == 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A'

        # get empty lines (where start and end the header)
        emptyLineIndex = []
        for lineIndex in range(len(lines)):
            if len(lines[lineIndex]) == 1:
                emptyLineIndex.append(lineIndex)

        # service information of the experiment, it should be in between the first two empty lines
        serviceInformation = []
        for ind in range(emptyLineIndex[0]+1, emptyLineIndex[1]+4):
            serviceInformation.append(get_metadata_in_one_line(lines[ind]))



        ###########  D A T A   A B O U T    E X P E R I M E N T  #######

        # humanreadable ID = name of experiment is the folder name of the data file
        metadata_emodule['humanreadableID'] = folderName  

        # ID of this experiment
        emoduleID = str(uuid.uuid4())
        metadata_emodule['ID'] = emoduleID

        # get experiment date and time in Protegé format YYYY-MM-DDTHH:mm:SS
        date = serviceInformation[10][4]  #datetime.datetime.strptime(,'%d.%m.%y')
        date_only = datetime.datetime.strptime(date.split(" ")[0], '%d.%m.%Y')
        date_protegeformat = date_only.strftime('%Y-%m-%d') + "T" + date.split(" ")[1]
        metadata_emodule['ExperimentDate'] = str(date_protegeformat)

        # operator name - This data has no placeholder yet.
        #metadata_emodule['tester_name'] = serviceInformation[2][1]

        # remarks - This data has no placeholder yet.
        #metadata_emodule['remark'] = serviceInformation[4][1]

        # set experiment lab location to BAM
        metadata_emodule['Lab'] = 'BAM'

        # set Compression and Transducer Column
        metadata_emodule['CompressionColumn'] = 0
        metadata_emodule['CompressionForce_Unit'] = "kN"
        metadata_emodule['TransducerColumn'] = [1, 2, 3]
        metadata_emodule['Extensometer_Unit'] = "mm"  # Transducer messen eine Verschiebung.
        metadata_emodule['CompressiveStrength'] = 55.8

        # set extensometer gauge length
        metadata_emodule['ExtensometerLength'] = 100
        metadata_emodule['ExtensometerLength_Unit'] = "mm"

        metadata_emodule['EModule_Unit'] = "GPa"  # laut Norm, mit einer Nachkommastelle; oft auch MPa oder N/mm²

        # set paths
        metadata_emodule['ProcessedFile'] = os.path.join('../usecases/MinimumWorkingExample/emodul/processed_data')  # path to csv file with values extracted by emodul_generate_processed_data.py
        metadata_emodule['RawDataFile'] = os.path.join(rawDataPath, specimen_file).replace('\\', '/')
        metadata_emodule['EModule'] = 33.06

  # path to specimen.dat
        try:
            #with open(str(rawDataPath), encoding="utf8", errors='ignore') as mix_data:
            with open(str(rawDataPath)+'/' + str(mix_file), encoding="utf8", errors='ignore') as mix_data:
            #with open(path_to_json) as mix_data:
                lines = mix_data.readlines()
                lines = lines[0].strip()

                #dataPath = Path(rawDataPath).parents[1]
                #metadata_emodule['MixDataFile']= os.path.join(dataPath, "Mischungen", lines)
        except:
            #metadata_emodule['MixDataFile'] = None
            raise Exception("No mixdesign json-file found!")


        ###########  D A T A   A B O U T    S P E C I M E N #######

        # name of specimen (humanreadable)
        #metadata_specimen['humanreadableID'] = serviceInformation[3][1]
        metadata_specimen['humanreadableID'] = folderName

        # ID of this specimen, save to specimen metadata and to emodule metadata
        #specimenID = str(uuid.uuid4())
        metadata_emodule['SpecimenID'] = metadata_specimen['ID'] = emoduleID

        # save Mixdesign ID to specimen metadata
        try:
            with open("../usecases/MinimumWorkingExample/mixture/metadata_json_files/" + lines[:-4] + "json", "r") as mixjson:  #change location of where
            #with open(path_to_json, "r") as mixjson:
                mixdesign = json.load(mixjson)
                #mixtureHumID = mixdesign['humanreadableID']
                #metadata_specimen['MixtureHumID'] = mixtureHumID
                mixtureID = mixdesign['ID']
                metadata_specimen['MixtureID'] = mixtureID
        except AttributeError:
            raise Exception("No mixdesign json-file found! Can't import the ID and save it to the output!")

        # set specimen age to 28 days
        metadata_emodule['SpecimenAge'] = 28.0
        metadata_emodule['SpecimenAge_Unit'] = 'day'

        # weight 
        metadata_specimen['SpecimenMass'] = float(replace_comma(serviceInformation[5][1]))
        metadata_specimen['SpecimenMass_Unit'] = 'g'

        # set size of specimen
        metadata_specimen['SpecimenDiameter'] = float(replace_comma(serviceInformation[6][1]))  #diameter
        metadata_specimen['SpecimenDiameter_Unit'] = 'mm'
        metadata_specimen['SpecimenLength'] = float(replace_comma(serviceInformation[7][1]))  #length
        metadata_specimen['SpecimenLength_Unit'] = 'mm'
        if metadata_specimen['SpecimenDiameter'] > metadata_specimen['SpecimenLength']:
            dir_name = metadata_emodule['ExperimentName']
            raise Exception(f'Diameter is larger than length, please fix the mistake in {dir_name}')


    return metadata_emodule, metadata_specimen


def emodul_metadata(rawDataPath, metaDataFile, specimenDataFile):
    """Creates two json files with extracted metadata, one for emodule and one
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
    
    with open(metaDataFile, 'w') as jsonFile:
        json.dump(metadata, jsonFile, sort_keys=False, ensure_ascii=False, indent=4)
    with open(specimenDataFile, 'w') as jsonFile:
        json.dump(specimen, jsonFile, sort_keys=False, ensure_ascii=False, indent=4)


def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from BAM e-module experiment')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata json
    parser.add_argument('-o', '--output', help='Path to extracted json files, one for emodule and one for specimen')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = '../usecases/MinimumWorkingExample/Data/E-modul/Wolf 8.2 Probe 1'
    if args.output == None:
        args.output = ['../usecases/MinimumWorkingExample/emodul/metadata_json_files/testMetaData.json', '../usecases/MinimumWorkingExample/emodul/metadata_json_files/testSpecimenData.json']



    # run extraction and write metadata file
    emodul_metadata(args.input, args.output[0], args.output[1])



if __name__ == "__main__":
    main()


