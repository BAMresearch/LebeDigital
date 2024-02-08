# Script for metadata-extraction for mixes with CEM I and CEM II. Output json.
# Works with the MixDesign ontology for mapping.


# ------------------------------------------------------------------------------

from cmath import nan
import pandas as pd
import numpy as np
# removing 'SettingWithCopyWarning: A value is trying to be set on a copy of a slice from a DataFrame'
pd.options.mode.chained_assignment = None  # default='warn'
import os
import json
from loguru import logger
from pathlib import Path
import argparse
import uuid

# Set up logger
baseDir = Path(__file__).parents[0]
logPath = os.path.join(baseDir, "logs", "file_{time}.log")
# logger.add(logPath, level="DEBUG")  # this also displays the log in the console
logger.configure(handlers=[{"sink": logPath, "level": "DEBUG"}])


# function to convert german formatting to english
def replace_comma(string, format='float'):
    if '---' in string:
        string = np.nan  # maybe None? But this will cause errors when float(string)
        return string
    elif format == 'float':
        string = string.replace(',', '.')
        return float(string)
    else:
        string = string.replace(',', '.')
        return string


# function to check for nan-values independently of the format (str/float)
def isNaN(num):
    return num is None or num != num


# decorater in case you want to catch errors so that the script won't break
# but just pass without output:
# @logger.catch

# extraction script
def extract_metadata_mixdesign(locationOfRawData):
    """
        Extracts the metadata from all "Rezeptur"-sheets of a given datafile
        (xls or xlsx). Creates also an ID. Returns a dictionary.

        Parameter
        ---------
        locationOfRawData : string
            Path of the excelsheet (xls or xlsx) containing the metadata in one
            or multiple "Rezeptur"-Sheet(s).

        Output
        -------
        The dict containing the metadata will be returned.
    """

    # Find sheets in the file containing the mixture (keyword: "Rezeptur"), allow
    # only one sheet per file
    excelsheet = os.path.basename(locationOfRawData)
    excelfile = pd.read_excel(locationOfRawData, sheet_name=None)
    listofkeys = [i for i in excelfile.keys() if 'Rezeptur' in i]
    logger.debug('Working on file: ' + excelsheet)
    logger.debug('Following sheet(s) contain mixture metadata in this file: ' + str(listofkeys))

    if len(listofkeys) != 1:
        logger.error('None or multiple sheets with mixture found in the raw data.')
        raise Exception('None or multiple sheets with mixture found in the raw data.')
    else:
        sheet = listofkeys[0]

        # name of json-file will be experiment-name
        name = os.path.basename(excelsheet).split('.xl')[0]

        # save data from excelsheet into pandas dataframe
        exceltodf = excelfile[sheet]

        # create empty dictionary for metadata
        metadata = {}

        # the layout of the Excel table can vary, the indices of labels are not
        # always the same; that's why: find now the indices of the labels and
        # store it in a dictionary
        labelidx = {}
        labelcolumn = exceltodf.iloc[:, 0]  # select first column (containing labels)
        for i in range(len(labelcolumn)):
            labelcolumn[i] = str(labelcolumn[i]).strip()  # remove whitespace

            # fill dictionary with labels and corresponding indices, unless the
            # label is "addition". Then differ between 1st and 2nd addition
            if labelcolumn[i] != 'Zusatzstoff':
                labelidx[labelcolumn[i]] = i
            elif labelcolumn[i] == 'Zusatzstoff' and 'Zusatzstoff1' not in labelidx.keys():
                labelidx['Zusatzstoff1'] = i
            elif labelcolumn[i] == 'Zusatzstoff' and 'Zusatzstoff1' in labelidx.keys() \
                    and 'Zusatzstoff2' not in labelidx.keys():
                labelidx['Zusatzstoff2'] = i
                logger.debug('Second addition found in raw data.')
            else:
                logger.error('More than two additions found in raw data.')
                raise Exception('More than two additions found in raw data.')

        # Check for missing labels; the following labels should exist (except
        # Zusatzstoff 2, not all raw files have two additions/Zusatzstoffe)
        default_labels = ['Bezeichnung der Proben:', 'Zement', 'Wasser (gesamt)',
                          'Zusatzmittel', 'Zuschlag (gesamt)', 'Zusatzstoff1', 'Zusatzstoff2']
        missing_labels = [i for i in default_labels if i not in labelidx.keys()]
        if len(missing_labels) != 0:
            if missing_labels == ['Zusatzstoff2']:
                logger.warning('No addition2 in raw data.')
            else:
                logger.error('Check raw data, there are labels missing: ' + str(missing_labels))
                raise KeyError('Check raw data, there are labels missing', missing_labels)

        # Some files don't have the type of addition/Zusatzstoff only labeled
        # in a cell that will be neglected during the extraction, so this saves
        # the type of Addition inside the annotation - but only in case it isn't
        # mentioned there already
        addition_finder = [True if i == 'Zusatzstoff' else False for i in labelcolumn]
        idx_addition = [i for i in range(len(addition_finder)) if addition_finder[i] == True]
        logger.debug('Number of additions in raw data: ' + str(len(idx_addition)))
        for i in idx_addition:
            # add the name in the annotation if not written there already
            if str(exceltodf.iloc[i, 1]) in str(exceltodf.iloc[i, 8]):
                pass
            elif isNaN(exceltodf.iloc[i, 8]):
                exceltodf.iloc[i, 8] = str(exceltodf.iloc[i, 1])
            else:
                exceltodf.iloc[i, 8] = str(exceltodf.iloc[i, 8]) + ' ' + str(exceltodf.iloc[i, 1])

        # This function will ensure that no empty annotation-information will
        # be passed to the json-file (check annotation-cell for nan)
        def no_empty_annotation(name):
            if isNaN(exceltodf.iat[idx, 8]):
                logger.debug('Empty annotation in ' + str(name))
                pass
            else:
                dic_label = str(name + '_Type')
                metadata[dic_label] = replace_comma(str(exceltodf.iat[idx, 8]), format='str')

        ############### E X T R A C T I O N #############

        # get raw data file name
        metadata['RawDataFile'] = locationOfRawData

        # get date (always the same position) & set time to 12:00 - Protege datetime format YYYY-MM-DDTHH:mm:SS
        metadata['MixingDate'] = str(exceltodf.columns[9])[:10] + "T12:00:00"

        # lab location - hardcoded
        metadata['Lab'] = "BAM"

        # ID of this mix
        MixID = str(uuid.uuid4())
        metadata['ID'] = MixID

        # humanreadable ID of this mix - is the name of the file without format
        metadata['humanreadableID'] = name

        # ----------------------------------------------------------------------

        # Extraction of the columns 'Stoffmenge' (QuantityInMix), 'Dichte bzw.
        # Rohdichte' (Density).

        # Cement data ('Zement')
        if 'Zement' not in missing_labels:
            idx = labelidx['Zement']
            metadata['Cement1_Content'] = float(replace_comma(str(exceltodf.iat[idx, 2])))
            metadata['Cement1_Content_Unit'] = 'kg/m^3'
            metadata['Cement1_Density'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Cement1_Density_Unit'] = 'kg/dm^3'
            no_empty_annotation('Cement1')
        else:
            logger.error('cement not included in json-file')

        # total water data ('Wasser (gesamt)')
        if 'Wasser (gesamt)' not in missing_labels:
            idx = labelidx['Wasser (gesamt)']
            metadata['Water_Content'] = float(replace_comma(str(exceltodf.iat[idx, 2])))
            metadata['Water_Content_Unit'] = 'kg/m^3'
            metadata['Water_Density'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Water_Density_Unit'] = 'kg/dm^3'
            no_empty_annotation('Water')
        else:
            logger.error('Water not included in json-file')

        # water cement ratio ('Wasserzementwert')
        try:
            metadata['WaterCementRatio'] = float(metadata['Water_Content']
                                                 / metadata['Cement1_Content'])
        except:
            raise Exception("Can not calculate water-cement-ratio! No values found!")

        # Admixture/Plasticizer ('Zusatzmittel')
        if 'Zusatzmittel' not in missing_labels:
            idx = labelidx['Zusatzmittel']
            metadata['Admixture1_Content'] = float(replace_comma(str(exceltodf.iat[idx, 2])))
            metadata['Admixture1_Content_Unit'] = 'kg/m^3'
            metadata['Admixture1_Density'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Admixture1_Density_Unit'] = 'kg/dm^3'
            no_empty_annotation('Admixture1')
        else:
            logger.error('Plasticizer/Admixture not included in json-file')

        # Aggregate ('Zuschlag (gesamt)')
        if 'Zuschlag (gesamt)' not in missing_labels:
            idx = labelidx['Zuschlag (gesamt)']
            metadata['Aggregate1_Content'] = float(replace_comma(str(exceltodf.iat[idx, 2])))
            metadata['Aggregate1_Content_Unit'] = 'kg/m^3'
            metadata['Aggregate1_Size'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Aggregate1_Size_Unit'] = float('nan')
            metadata['Aggregate1_Density'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Aggregate1_Density_Unit'] = 'kg/dm^3'
            no_empty_annotation('Aggregate1')
        else:
            logger.error('Okrilla/aggregate not included in json-file')

        # Addition data ('Zusatzstoff')
        if 'Zusatzstoff1' not in missing_labels:
            idx = labelidx['Zusatzstoff1']
            metadata['Addition1_Content'] = float(replace_comma(str(exceltodf.iat[idx, 2])))
            metadata['Addition1_Content_Unit'] = 'kg/m^3'
            metadata['Addition1_Density'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Addition1_Density_Unit'] = 'kg/dm^3'
            no_empty_annotation('Addition1')
        else:
            logger.error('addition not included in json-file')

        # Addition data ('Zusatzstoff')
        if 'Zusatzstoff2' not in missing_labels:
            idx = labelidx['Zusatzstoff2']
            metadata['Addition2_Content'] = float(replace_comma(str(exceltodf.iat[idx, 2])))
            metadata['Addition2_Content_Unit'] = 'kg/m^3'
            metadata['Addition2_Density'] = float(replace_comma(str(exceltodf.iat[idx, 4])))
            metadata['Addition2_Density_Unit'] = 'kg/dm^3'
            no_empty_annotation('Addition2')
        else:
            logger.error('addition2 not included in json-file')

        return metadata


def mix_metadata(rawDataPath, metaDataFile):
    """Creates a json file with extracted metadata for the mixDesign.

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    metaDataFile : string
        Path to the output data file for mix metadata
    """

    # extracting the metadata
    metadata = extract_metadata_mixdesign(rawDataPath)
    # Replace occurrences of NaN with None in the metadata dictionary
    metadata = {key: None if isNaN(value) else value for key, value in metadata.items()}

    json_name = rawDataPath.name.split('.')[0]
    metaDataFile = str(metaDataFile) + json_name
    # print(rawDataPath.split('/')[-1].split('.')[0])
    # writing the metadata to json file
    with open(metaDataFile + ".json", 'w') as jsonFile:
        json.dump(metadata, jsonFile, sort_keys=False, ensure_ascii=False, indent=4)

    return metaDataFile


def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from MixDesign.')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata json
    parser.add_argument('-o', '--output', help='Path to extracted json files.')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = '../../usecases/MinimumWorkingExample/Data/Mischungen/2019_06_26 Klimek Geschossdecke_Quarzkies.xls'
        
    if args.output == None:
        args.output = '../../usecases/MinimumWorkingExample/mixture/metadata_json_files/'
       
    # run extraction and write metadata file
    # path_to_json = mix_metadata(args.input, args.output)
    mix_metadata(args.input, args.output)

    # return path_to_json


if __name__ == "__main__":
    main()
