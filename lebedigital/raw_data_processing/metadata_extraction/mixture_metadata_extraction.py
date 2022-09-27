# Script for metadata-extraction: 
# give path to raw data and output-path, get yaml-file with the metadata
# If the raw data (excelfile) has multiple sheets, all of them are processed.
# based on extraction-code of Erik


#                 !!! Current Issues / To-Do !!!
# The script only works on files that have exactly the same structure/labels as 
# "2014_08_05 Rezeptur_MI.xlsx". Information about the kind of Zusatzstoff is 
# missing, also multiple Zusatzstoffe can't be processed. 
# Also I am not sure about all translations, I am waiting for an update of the 
# taxonomy-file ("Luftgehalt", "Zuschlag",....).


# Workflow:
# 1. Installing and importing the necessary libraries
# 2. Define a function to convert german formatting to english
# 3. Defining a function to extract metadata:
# 3.a Look for all sheets containing the word "Rezeptur" and load them
# 3.b For each sheet: find indices of the labels
# 3.c For each sheet: create and fill a dic & export dic as yaml-file

#------------------------------------------------------------------------------

from cmath import nan
import pandas as pd
import os
import yaml
from loguru import logger 

# function to convert german formatting to english
def replace_comma(string):
    if "---" in string:
        string = nan  # maybe None? But this will cause errors when float(string)
        return string
    else:
        string = string.replace(',', '.')
        return float(string)


# extraction script
def extract_metadata_mixture(
        locationOfRawData,
        locationOfProcessedData
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
            Path of the target folder for yaml-file.

        Output
        -------
        yamlFile containing the dictionary with the metadata
                    
    """
    

    # Find sheets in the file containing the mixture (keyword: "Rezeptur")
    excelsheet = os.path.basename(locationOfRawData)
    logger.debug("Working on: "+ excelsheet)
    excelfile = pd.read_excel(locationOfRawData, sheet_name= None) 
    listofkeys = [i for i in excelfile.keys() if "Rezeptur" in i] 
    logger.debug("Sheets containing mixture metadata in this file: " + str(listofkeys))


    for sheet in listofkeys:

        ##############  S E T U P ##############

        # name of yaml-file will be experiment-name + sheet name
        name = os.path.basename(excelsheet).split(".xl")[0] + " ___ " + sheet
        
        # save data from excelsheet into pandas dataframe
        exceltodf = excelfile[sheet]
     
        # create empty dictionary for metadata
        metadata = {}

        # the layout of the excel table can vary, the indices of labels are not 
        # always the same; that's why: find now the indices of the labels
        labelidx = {}
        labelcolumn = exceltodf.iloc[:,0]  # select first column (containing labels)
        logger.debug("All information of the label column:" + labelcolumn[0])
        for i in range(len(labelcolumn)):
            labelcolumn[i] = str(labelcolumn[i]).strip()  # remove whitespace
            labelidx[labelcolumn[i]] = i  # fill dictionary
        logger.debug(labelidx)

        # check for missing and additional labels; the following labels have to exist
        default_labels = ['Bezeichnung der Proben:', 'Zement', 'Wasser (gesamt)', 'Wasser (wirksam)', 'Luftgehalt', 
                        'Zusatzstoff', 'Zusatzmittel', 'Zuschlag (gesamt)',
                        'Frischbeton', 'Mehlkornanteil', 'Mörtelanteil']
        missing_labels =  [i for i in default_labels if i not in labelidx.keys()]
        logger.debug("Labels missing in the raw data: ",  missing_labels)
        # these are all labels of "2014_08_05 Rezeptur_MI.xlsx"
        # ['7.0', 'nan', 'Antragsteller:', 'Antrags-/ Projekt-Nr.:', 'Bezeichnung der Proben:', 
        # 'Betonsorte u.-festigkeitsklasse:', 'Wasserzementwert', 'Sieblinie n. DIN 1045:', 
        # 'Sieblinie des Zuschlags:', 'Sieblochweite in mm', 'Durchgang in Vol. -%', 
        # 'Berechnung der Betonzusammensetzung', 'Stoffart', 'Zement', 'Wasser (gesamt)', 
        # 'Wasser (wirksam)', 'Luftgehalt', 'Zusatzstoff', 'Zusatzmittel', 'Gesamt', 'Zuschlag (gesamt)', 
        # '0 / 0,3', '0,1 / 0,5', '0,5 / 1,0', '1,0 / 2,0', '2,0 / 4,0', '4,0 / 8,0', '8,0 / 16,0', 
        #'Frischbeton', 'Mehlkornanteil', 'Mörtelanteil']
        

        # Some files have multiple entries of Zusatzstoffe, which are not further labeled, so this adds 
        # the type of Zusatzstoff
        howmanyzusatzstoffe = [True if i == word else False for i in exceltodf["Stoffart"] ]
        indexzusatzstoff = [i for i in range(len(howmanyzusatzstoffe)) if howmanyzusatzstoffe[i] == True]
        for i in indexzusatzstoff:
            exceltodf.iloc[i,0] += " " + str(exceltodf.iloc[i,1])



        ############### E X T R A C T I O N #############

        # get date and time (always the same position)
        metadata['operator_date'] = str(exceltodf.columns[9])[:10]

        # operator name (always the same position)
        metadata['tester_name'] = exceltodf.iloc[0,10]

        # name of specimen
        idx = labelidx["Bezeichnung der Proben:"]
        metadata['specimen_name'] = exceltodf.iloc[idx,3] #4

        #----------------------------------------------------------------------

        # Extraction of the columns "Stoffmenge" (QuantityInMix), "Dichte bzw. 
        # Rohdichte" (BulkDensity), "Stoffraum" (Volume), "Sonstiges / Bemerkungen"
        # (Annotation):

        # Cement data ("Zement") # idx 20 for 2014_08_05 Rezeptur_MI
        idx = labelidx["Zement"]
        metadata['cement--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
        metadata['cement--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        metadata['cement--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
        metadata['cement--Annotation'] = str(exceltodf.iat[idx,8])

        # total water data ("Wasser (gesamt)") # idx 21 for 2014_08_05 Rezeptur_MI
        idx = labelidx["Wasser (gesamt)"]
        metadata['water_total--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
        metadata['water_total--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        metadata['water_total--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
        metadata['water_total--Annotation'] = str(exceltodf.iat[idx,8])

        # water cement ratio ("Wasserzementwert")
        metadata["watercementratio"] = float(metadata['water_total--QuantityInMix'] / metadata['cement--QuantityInMix'])

        # effective water data ("Wasser (wirksam)") # idx 22 for 2014_08_05 Rezeptur_MI 
        idx = labelidx["Wasser (wirksam)"]
        metadata['water_effective--QuantityInMix'] = replace_comma(str(exceltodf.iat[idx,2]))
        metadata['water_effective--BulkDensity'] = replace_comma(str(exceltodf.iat[idx,4]))
        metadata['water_effective--Volume'] = replace_comma(str(exceltodf.iat[idx,6]))
        metadata['water_effective--Annotation'] = str(exceltodf.iat[idx,8])

        # air content data ("Luftgehalt") # idx 23 for 2014_08_05 Rezeptur_MI
        idx = labelidx["Luftgehalt"]
        metadata['air_content--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
        metadata['air_content--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        metadata['air_content--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
        metadata['air_content--Annotation'] = str(exceltodf.iat[idx,8])

        # Admixture data ("Zusatzstoff") # idx 24 for 2014_08_05 Rezeptur_MI
        idx = labelidx["Zusatzstoff"]
        metadata['admixture--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
        metadata['admixture--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        metadata['admixture--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
        metadata['admixture--Annotation'] = str(exceltodf.iat[idx,8])

        # Zusatzmittel MISSING
        idx = labelidx["Zusatzmittel"]

        # Aggregate data ("Zuschlag (gesamt)") # idx 27 for 2014_08_05 Rezeptur_MI
        idx = labelidx["Zuschlag (gesamt)"]
        metadata['zuschlag--QuantityInMix'] = float(replace_comma(str(exceltodf.iat[idx,2])))
        metadata['zuschlag--BulkDensity'] = float(replace_comma(str(exceltodf.iat[idx,4])))
        metadata['zuschlag--Volume'] = float(replace_comma(str(exceltodf.iat[idx,6])))
        metadata['zuschlag--Annotation'] = str(exceltodf.iat[idx,8])

        # ("Frischbeton") MISSING
        #idx = labelidx["Frischbeton"]

        # ("Mehlkornanteil") MISSING
        #idx = labelidx["Mehlkornanteil"]

        # ("Mörtelanteil") MISSING
        #idx = labelidx["Mörtelanteil"]

        #print(metadata)



        ############### O U T P U T #############
 
        # # Write the data to a yaml file
        with open(os.path.join(locationOfProcessedData, name + '.yaml'), mode='w') as yamlFile:
           yaml.dump(metadata, yamlFile, sort_keys=False)
           

extract_metadata_mixture("C:\\Users\\afriemel\\Documents\\Backup\\mixture\\2014_08_05 Rezeptur_MI.xlsx","C:\\Users\\afriemel\\Documents\\Backup\\mixture")