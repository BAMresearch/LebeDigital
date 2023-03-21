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
def extract_metadata_mixdesign(
        locationOfRawData,
        locationOfProcessedData = None
        ):    

    """
        Extracts the metadata from all "Rezeptur"-sheets of a given datafile 
        (xls or xlsx). Creates one yaml-file per sheet containing the keyword
        "Rezeptur".


        Parameter
        ---------
        locationOfRawData : string
            Path of the excelsheet (xls or xlsx) containing the metadata in one 
            or multiple "Rezeptur"-Sheet(s).
        locationOfProcessedData : string
            Path of the target folder for yaml-file (optional, give only if you 
            want a yaml-file to be generated).

        Output
        -------
        If no ouput path (locationOfProcessedData) is given (f.e. for unittesting), 
        the dict containing the metadata will be returned. Otherwise a yaml file 
        will be created.
                    
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

        # get date and time (always the same position)
        metadata['MixingDate'] = str(exceltodf.columns[9])[:10]

        # lab location - hardcoded
        metadata['Lab'] = "BAM"

        #----------------------------------------------------------------------

        # Extraction of the columns 'Stoffmenge' (QuantityInMix), 'Dichte bzw. 
        # Rohdichte' (Density).

        # Cement data ('Zement') 
        if 'Zement' not in missing_labels:
            idx = labelidx['Zement']
            metadata['CEMIQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['CEMIDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        else:
            logger.error('cement not included in yaml-file')


        # total water data ('Wasser (gesamt)') 
        if 'Wasser (gesamt)' not in missing_labels:
            idx = labelidx['Wasser (gesamt)']
            metadata['MixingWaterQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['WaterDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
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
            metadata['AirDensity'] = float(0)
        else:
            logger.error('AirContent not included in yaml-file')


        # Admixture/Plasticizer ('Zusatzmittel') 
        if 'Zusatzmittel' not in missing_labels:
            idx = labelidx['Zusatzmittel']
            metadata['PlasticizerQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['PlasticizerDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        else:
            logger.error('Plasticizer/Admixture not included in yaml-file')


        # Aggregate ('Zuschlag (gesamt)')
        if 'Zuschlag (gesamt)' not in missing_labels:
            idx = labelidx['Zuschlag (gesamt)']
            metadata['OkrillaQtyInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['OkrillaDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        else:
            logger.error('Okrilla/aggregate not included in yaml-file')



        ############################ O U T P U T #############################
        if locationOfProcessedData == None:
            return metadata
        else:
            with open(os.path.join(locationOfProcessedData, name + '.yaml'), mode='w') as yamlFile:
                yaml.dump(metadata, yamlFile, sort_keys=False, allow_unicode=True)
           

