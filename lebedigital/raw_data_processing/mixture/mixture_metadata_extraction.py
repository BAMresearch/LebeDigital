# Script for metadata-extraction: 
# give path to raw data and output-path, get yaml-file with the metadata
# If the raw data (excelfile) has multiple sheets, all of them are processed.
# based on extraction-code of Erik


# Workflow:
# 1. Installing and importing the necessary libraries
# 2. Define a function to convert german formatting to english
# 3. Define a function to check for nan value regardless of float or string
# 4. Defining a function to extract metadata:
# 3.a Look for all sheets containing the word "Rezeptur" in the file, load them
#     and go through them one by one
# 3.b Find indices of the labels, except for labels called "addition" (Zusatzstoff)
#     (those get renamed to addition1 and addition2 to ensure unique keys), and
#     create a dictionary with this information
# 3.c Check for missing labels
# 3.d Save the name of the addition into the annotation in case it isn't mentioned
#     there already
# 3.e Define a function to ignore empty annotations.
# 3.f Extract the metadata for each label which is existing. 
# 3.g Depending on wheter an output-path is given or not, return the metadata
#     dictionary or create a yaml-file.

#------------------------------------------------------------------------------

from cmath import nan
import pandas as pd
import os
import yaml
from loguru import logger 
from pathlib import Path


# Set up logger
baseDir = Path(__file__).parents[0]
logPath = os.path.join(baseDir, "logs","file_{time}.log")
logger.add(logPath, level="DEBUG")


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
def extract_metadata_mixture(
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
    

    # Find sheets in the file containing the mixture (keyword: "Rezeptur")
    excelsheet = os.path.basename(locationOfRawData)
    excelfile = pd.read_excel(locationOfRawData, sheet_name= None) 
    listofkeys = [i for i in excelfile.keys() if 'Rezeptur' in i] 
    logger.debug('Working on file: '+ excelsheet)
    logger.debug('Following sheets contain mixture metadata in this file: ' + str(listofkeys))


    for sheet in listofkeys:

        logger.debug('Working on sheet: '+ sheet)
        ##############  S E T U P ##############

        # name of yaml-file will be experiment-name + sheet name
        name = os.path.basename(excelsheet).split('.xl')[0] + ' ___ ' + sheet
        
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
                        'Wasser (wirksam)', 'Luftgehalt', 'Zusatzstoff1', 'Zusatzstoff2',
                        'Zusatzmittel', 'Zuschlag (gesamt)'] 
        missing_labels =  [i for i in default_labels if i not in labelidx.keys()]
        if len(missing_labels) != 0:
            if missing_labels == ['Zusatzstoff2']:
                logger.warning('No addition2 in raw data.')
            else:
                logger.error('Check raw data, there are labels missing.')
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
            if str(exceltodf.iloc[i,1]) in str(exceltodf.iloc[i,8]):
                pass
            elif isNaN(exceltodf.iloc[i,8]):
                exceltodf.iloc[i,8] =  str(exceltodf.iloc[i,1])
            else:
                exceltodf.iloc[i,8] =  str(exceltodf.iloc[i,8]) +' - ' + str(exceltodf.iloc[i,1])


        # This function will ensure that no empty annotation-information will 
        # be passed to the yaml-file (check annotation-cell for nan)
        def no_empty_annotation(name):
            if isNaN(exceltodf.iat[idx,8]):
                logger.debug('Empty annotation in ' + str(name))
                pass
            else:
                dic_label = str(name + '--Annotation')
                metadata[dic_label] = replace_comma(str(exceltodf.iat[idx,8]), format='str')


        ############### E X T R A C T I O N #############

        # get date and time (always the same position)
        metadata['operator_date'] = str(exceltodf.columns[9])[:10]

        # operator name (always the same position)
        metadata['tester_name'] = exceltodf.iloc[0,10]

        # name of specimen
        idx = labelidx['Bezeichnung der Proben:']
        metadata['specimen_name'] = exceltodf.iloc[idx,3]

        #----------------------------------------------------------------------

        # Extraction of the columns 'Stoffmenge' (QuantityInMix), 'Dichte bzw. 
        # Rohdichte' (BulkDensity), 'Stoffraum' (Volume), 'Sonstiges / Bemerkungen'
        # (Annotation):

        # Cement data ('Zement') 
        if 'Zement' not in missing_labels:
            idx = labelidx['Zement']
            metadata['cement--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['cement--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['cement--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('cement')
        else:
            logger.error('cement not included in yaml-file')

        # total water data ('Wasser (gesamt)') 
        if 'Wasser (gesamt)' not in missing_labels:
            idx = labelidx['Wasser (gesamt)']
            metadata['water_total--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['water_total--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['water_total--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('water_total')
        else:
            logger.error('water_total not included in yaml-file')


        # water cement ratio ('Wasserzementwert')
        if 'Zement' not in missing_labels and 'Wasser (gesamt)' not in missing_labels:
            metadata['water_cement_ratio'] = float(metadata['water_total--QuantityInMix'] 
                                                    / metadata['cement--QuantityInMix'])
        else:
            logger.error('water_cement_ratio not included in yaml-file')


        # effective water data ('Wasser (wirksam)')  
        if 'Wasser (wirksam)' not in missing_labels:
            idx = labelidx['Wasser (wirksam)']
            metadata['water_effective--QuantityInMix'] = replace_comma(str(exceltodf.iat[idx,2]))
            metadata['water_effective--BulkDensity'] = replace_comma(str(exceltodf.iat[idx,4]))
            metadata['water_effective--Volume'] = replace_comma(str(exceltodf.iat[idx,6]))
            no_empty_annotation('water_effective')
        else:
            logger.error('water_effective not included in yaml-file')


        # air content data ('Luftgehalt') 
        if 'Luftgehalt' not in missing_labels:
            idx = labelidx['Luftgehalt']
            metadata['air_content--QuantityInMix'] = float(0)
            metadata['air_content--BulkDensity'] = float(0)
            metadata['air_content--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('air_content')
        else:
            logger.error('air_content not included in yaml-file')


        # Addition data ('Zusatzstoff') 1 
        if 'Zusatzstoff1' not in missing_labels:
            idx = labelidx['Zusatzstoff1']
            metadata['addition1--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['addition1--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['addition1--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('addition1')
        else:
            logger.error('addition1 not included in yaml-file')


        # Addition data ('Zusatzstoff') 2 
        if 'Zusatzstoff2' not in missing_labels:
            idx = labelidx['Zusatzstoff2']
            metadata['addition2--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['addition2--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['addition2--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('addition2')
            print( metadata['addition2--QuantityInMix'] )
        else:
            logger.warning('addition2 not included in yaml-file')


        # Admixture ('Zusatzmittel') 
        if 'Zusatzmittel' not in missing_labels:
            idx = labelidx['Zusatzmittel']
            metadata['admixture--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['admixture--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['admixture--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('admixture')
        else:
            logger.error('admixture not included in yaml-file')


        # Aggregate ('Zuschlag (gesamt)')
        if 'Zuschlag (gesamt)' not in missing_labels:
            idx = labelidx['Zuschlag (gesamt)']
            metadata['aggregate--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
            metadata['aggregate--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
            metadata['aggregate--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
            no_empty_annotation('aggregate')
        else:
            logger.error('aggregate not included in yaml-file')



        ############################ O U T P U T #############################
        if locationOfProcessedData == None:
            return metadata
        else:
            with open(os.path.join(locationOfProcessedData, name + '.yaml'), mode='w') as yamlFile:
                yaml.dump(metadata, yamlFile, sort_keys=False, allow_unicode=True)
           

