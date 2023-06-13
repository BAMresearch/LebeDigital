# Script for metadata-extraction for mixes with CEM I and CEM II. Output yaml
# should work with the MixDesign ontology for mapping.


#------------------------------------------------------------------------------

from cmath import nan
import pandas as pd
# removing 'SettingWithCopyWarning: A value is trying to be set on a copy of a slice from a DataFrame'
pd.options.mode.chained_assignment = None  # default='warn'

import os
import yaml
from loguru import logger 
from pathlib import Path
import argparse
import uuid


# Set up logger
baseDir = Path(__file__).parents[0]
logPath = os.path.join(baseDir, "logs","file_{time}.log")
#logger.add(logPath, level="DEBUG")  # this also displays the log in the console
logger.configure(handlers=[{"sink": logPath, "level": "DEBUG"}])



# function to convert german formatting to english
def replace_comma(string, format = 'float'):
    if '---' in string:
        string = nan  # maybe None? But this will cause errors when float(string)
        return string
    elif format == 'float':
        string = string.replace(',', '.')
        return float(string)
    else:
        string = string.replace(',', '.')
        return string



# function to check for nan-values independently of the format (str/float)
def isNaN(num):
    return num!= num



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
    excelfile = pd.read_excel(locationOfRawData, sheet_name= None) 
    listofkeys = [i for i in excelfile.keys() if 'Rezeptur' in i] 
    logger.debug('Working on file: '+ excelsheet)
    logger.debug('Following sheet(s) contain mixture metadata in this file: ' + str(listofkeys))

    if len(listofkeys) != 1:
        logger.error('None or multiple sheets with mixture found in the raw data.')
        raise Exception('None or multiple sheets with mixture found in the raw data.')
    else:
        sheet = listofkeys[0]

        # name of yaml-file will be experiment-name 
        name = os.path.basename(excelsheet).split('.xl')[0]
        
        # save data from excelsheet into pandas dataframe
        exceltodf = excelfile[sheet]
     
        # create empty dictionary for metadata
        metadata = {}

        # the layout of the excel table can vary, the indices of labels are not 
        # always the same; that's why: find now the indices of the labels and 
        # store it in a dictionary
        labelidx = {}
        labelcolumn = exceltodf.iloc[:,0]  # select first column (containing labels)
        for i in range(len(labelcolumn)):
            labelcolumn[i] = str(labelcolumn[i]).strip()  # remove whitespace
            labelidx[labelcolumn[i]] = i

        # Check for missing labels; the following labels should exist (except 
        # Zusatzstoff 2, not all raw files have two additions/Zusatzstoffe)
        default_labels = ['Bezeichnung der Proben:', 'Zement', 'Wasser (gesamt)', 
                        'Luftgehalt', 'Zusatzmittel', 'Zuschlag (gesamt)'] 
        missing_labels =  [i for i in default_labels if i not in labelidx.keys()]
        if len(missing_labels) != 0:
            logger.error('Check raw data, there are labels missing: ' + str(missing_labels))
            raise KeyError('Check raw data, there are labels missing', missing_labels)


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

        # human readable ID of this mix - is the name of the file without format
        metadata['humanreadableID'] = name

        #----------------------------------------------------------------------

        # Extraction of the columns 'Stoffmenge' (QuantityInMix), 'Dichte bzw. 
        # Rohdichte' (Density).

        # Cement data ('Zement') 
        if 'Zement' not in missing_labels:
            idx = labelidx['Zement']
            metadata['CEMIQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['CEMIQtyInMix_Unit'] = 'kg/m^3'
            metadata['CEMIDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['CEMIDensity_Unit'] = 'kg/dm^3'
        else:
            logger.error('cement not included in yaml-file')


        # total water data ('Wasser (gesamt)') 
        if 'Wasser (gesamt)' not in missing_labels:
            idx = labelidx['Wasser (gesamt)']
            metadata['MixingWaterQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['MixingWaterQtyInMix_Unit'] = 'kg/m^3'
            metadata['WaterDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['WaterDensity_Unit'] = 'kg/dm^3'
        else:
            logger.error('Water not included in yaml-file')


        # water cement ratio ('Wasserzementwert')
        if 'Zement' not in missing_labels and 'Wasser (gesamt)' not in missing_labels:
            metadata['WaterCementRatio'] = float(metadata['MixingWaterQtyInMix'] 
                                                    / metadata['CEMIQtyInMix'])
        else:
            logger.error('WaterCementRatio not included in yaml-file')


        # air content ('Luftgehalt') 
        if 'Luftgehalt' not in missing_labels:
            idx = labelidx['Luftgehalt']
            metadata['AirContent'] = float(0) # Quantity
            metadata['AirContent_Unit'] = 'kg/m^3'
            metadata['AirDensity'] = float(0)
            metadata['AirDensity_Unit'] = 'kg/dm^3'
        else:
            logger.error('AirContent not included in yaml-file')


        # Admixture/Plasticizer ('Zusatzmittel') 
        if 'Zusatzmittel' not in missing_labels:
            idx = labelidx['Zusatzmittel']
            metadata['PlasticizerQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['PlasticizerQtyInMix_Unit'] = 'kg/m^3'
            metadata['PlasticizerDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['PlasticizerDensity_Unit'] = 'kg/dm^3'
        else:
            logger.error('Plasticizer/Admixture not included in yaml-file')


        # Aggregate ('Zuschlag (gesamt)')
        if 'Zuschlag (gesamt)' not in missing_labels:
            idx = labelidx['Zuschlag (gesamt)']
            metadata['OkrillaQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['OkrillaQtyInMix_Unit'] = 'kg/m^3'
            metadata['OkrillaDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['OkrillaDensity_Unit'] = 'kg/dm^3'
        else:
            logger.error('Okrilla/aggregate not included in yaml-file')


        return metadata



def mix_metadata(rawDataPath, metaDataFile):
    """Creates a yaml file with extracted metadata for the mixDesign.

    Parameters
    ----------
    rawDataFile : string
        Path to the raw data file
    metaDataFile : string
        Path to the output data file for mix metadata
    """

    # extracting the metadata
    metadata = extract_metadata_mixdesign(rawDataPath)

    # writing the metadata to yaml file
    with open(metaDataFile, 'w') as yamlFile:
        yaml.dump(metadata, yamlFile, sort_keys=False, allow_unicode=True)


def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from MixDesign.')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata yaml
    parser.add_argument('-o', '--output', help='Path to extracted yaml files.')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = '../../../usecases/MinimumWorkingExample/Data/Mischungen/2014_08_05 Rezeptur_MI.xlsx'
    if args.output == None:
        args.output = '../../../usecases/MinimumWorkingExample/mixture/metadata_yaml_files/2014_08_05 Rezeptur_MI.yaml'

    # run extraction and write metadata file
    mix_metadata(args.input, args.output)

if __name__ == "__main__":
    main()



