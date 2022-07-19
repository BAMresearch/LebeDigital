# Script for metadata-extraction: function takes in as argument a given path to
# raw data and an output-path and returns a yaml-file with the metadata
#
# Workflow:
# 1. Installing and importing the necessary libraries
# 2. Defining a function to translate the labels from german to english
# 3. Defining a function to extract metadata:
# 3.a Look for sheets containing the word "Rezeptur" and load them
# 3.b Reduce the dataframe to only necessary entries
# 3.c Merge dataframe to a series
# 3.d Export series as yaml-file
#----------------------------------------------------------------------------------------------------------

# In[1]: Setup

from cmath import nan
import pandas as pd
import os
import yaml

def translate(
        df, 
        dic
        ):

    """
        Translates german labels to english according to the translation.xls

        Parameter
        ---------
        df
            The dataframe with the german labels that should be translated.
        dic
            The excelfile that contains the translations, columns are "german" and "english"
    """
    
    translation = pd.read_excel(dic)
    translation.columns = ["german", "english"]

    for deutsch, english in translation.itertuples(index=False,name=None):
        df = df.replace(deutsch,english,regex=True)
    return df



def extraction(
        locationOfRawData,
        locationOfProcessedData,
        english = True
        ):    

    """
        Extractes the metadata from the "Rezeptur"-Sheet of a given datafile (xls or xlsx).

        Parameter
        ---------
        locationOfRawData
            Path of the excelsheet (xls or xlsx) containing the metadata in one or multiple "Rezeptur"-Sheet(s).
        locationOfProcessedData
            Path of the target folder for yaml-file.
        english
            Boolean, checking if the metadata labels should be english or german 
            -> For True, it uses a provided dictionary. 
        
    """
    
    excelsheet = os.path.basename(locationOfRawData)
    print("\n Working on: "+ excelsheet)
    excelfile = pd.read_excel(locationOfRawData, sheet_name= None) 
    listofkeys = [i for i in excelfile.keys() if "Rezeptur" in i] 
    print("Relevant sheets in this file: " + str(listofkeys))


    for sheet in listofkeys:

        name = os.path.basename(excelsheet).split(".xl")[0] + " ___ " + sheet
        exceltodf = excelfile[sheet]
        
        # extract header, translate and set proper index & column headers
        exceltodf.iat[17,4] = "Dichte"
        if english == True:
            exceltodf = translate(exceltodf,"translation.xlsx")
            specimeninfo = {
                    "Date" : str(exceltodf.columns[9])[:10], 
                    "Editor" : exceltodf.iloc[0,10], 
                    "Requester" : exceltodf.iloc[2,3],
                    "Project number": exceltodf.iloc[3,3],
                    "Specimen" : exceltodf.iloc[4,3], 
                    "Betonsorte u -festigkeitsklasse" : exceltodf.iloc[6,4],
                    "Wasserzementwert" : exceltodf.iloc[7,2],
                    "Konsistenz" : exceltodf.iloc[7,7],
                    "Sieblinie n. DIN 1045" : exceltodf.iloc[8,2],
                    "Körnungsziffer" : exceltodf.iloc[8,8],
                    "Vol Leim/Zuschlag" : exceltodf.iloc[10,10]
                    }
            word = "Supplementary cementious materials"
        else:
            specimeninfo = {
                    "Datum" : str(exceltodf.columns[9])[:10], 
                    "Bearbeiter" : exceltodf.iloc[0,10], 
                    "Antragsteller" : exceltodf.iloc[2,3],
                    "Antrags-/Projektnummer": exceltodf.iloc[3,3],
                    "Bezeichnung der Proben" : exceltodf.iloc[4,3], 
                    "Betonsorte u -festigkeitsklasse" : exceltodf.iloc[6,4],
                    "Wasserzementwert" : exceltodf.iloc[7,2],
                    "Konsistenz" : exceltodf.iloc[7,7],
                    "Sieblinie n. DIN 1045" : exceltodf.iloc[8,2],
                    "Körnungsziffer" : exceltodf.iloc[8,8],
                    "Vol Leim/Zuschlag" : exceltodf.iloc[10,10]
                    }
            word = "Zusatzstoff"
        specimeninfo = pd.Series(specimeninfo)
        specimeninfo = specimeninfo.fillna("---") # fill nans
        exceltodf.rename(columns=exceltodf.iloc[17], inplace = True) # set column header
        exceltodf["Stoffart"] = exceltodf["Stoffart"].str.strip()  # remove whitespace
        exceltodf["Sonstiges/Bemerkungen"] = exceltodf["Sonstiges/Bemerkungen"].str.strip()  # remove whitespace
        
        # Some files have multiple entries of Zusatzstoffe, which are not further labeled, so this adds 
        # the type of Zusatzstoff
        howmanyzusatzstoffe = [True if i == word else False for i in exceltodf["Stoffart"] ]
        indexzusatzstoff = [i for i in range(len(howmanyzusatzstoffe)) if howmanyzusatzstoffe[i] == True]
        for i in indexzusatzstoff:
            exceltodf.iloc[i,0] += " " + str(exceltodf.iloc[i,1])
     
        # create new dataframe with only relevant data and chose Index column
        relevant_data = exceltodf.iloc[20:,[0,2,4,8]]
        relevant_data = relevant_data.set_index("Stoffart")
        
        # replace NaN 
        relevant_data[relevant_data.columns[-1]] = relevant_data[relevant_data.columns[-1]].fillna("---")
        relevant_data = relevant_data.drop(nan)

        # merge df into a series
        df = relevant_data
        df_out = df.stack()
        df_out.index = df_out.index.map('{0[0]} -- {0[1]}'.format)
        df_out = pd.concat([df_out, specimeninfo], ignore_index=False)

         # Write the data to a yaml file
        with open(os.path.join(locationOfProcessedData, name + '.yaml'), mode='w') as file:
            yaml.dump(df_out.to_dict(), file, sort_keys=False, allow_unicode=True, width=72, indent=4)

