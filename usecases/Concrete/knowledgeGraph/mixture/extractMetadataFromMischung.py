# Script for metadata-extraction
#
# Workflow:
# 1. Installing and importing the necessary libraries
# 2. Checking which excel-files exist and need to be extracted
# 3. Load a an excel file and check for the sheets
# 4. Convert each relevant sheet of each excelfile into one yaml-file

#----------------------------------------------------------------------------------------------------------

# In[1]: Setup

from cmath import nan
import pandas as pd
from IPython.display import display
import glob as g
import os
import yaml

# Locate the mixture data (concrete > Data > Mischungen)
datapath = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', 'Data', 'Mischungen'))
outputpath = os.path.realpath(os.path.join(os.path.dirname(__file__), 'Output'))

# Displaying setting
pd.set_option('display.max_rows', None)

#In[2]: Data collection
# Collect a list of excel-Files in the current folder
excelsheetlistall = g.glob(datapath + '/*xlsx')  + g.glob(datapath + './*xls') # Implement later, causes issues right now
excelsheetlist = [i for i in excelsheetlistall if "~" not in i] 
datanamelist = []
for file in excelsheetlist:
    datanamelist.append(os.path.basename(file))
    #datanamelist.append(os.path.basename(file).split(".xl")[0])
print("Found " + str(len(datanamelist)) + " excel-sheets in the directory: " + "\n"+ str(datanamelist) + "\n")

# Check which sheets have been already extracted
yamllist = g.glob(outputpath + './*yaml')
alreadydone = []
alreadydone.clear()
for file in yamllist:
    alreadydone.append(os.path.basename(file).split(".yaml")[0])

print("These " + str(len(alreadydone)) + " yaml-files already exist: " + "\n"+ str(alreadydone) + "\n")
todo = [i for i in datanamelist if any(i[:-5] in j for j in alreadydone)==False]#i.find(i) != -1 ]
print("These " + str(len(todo)) + " files will be processed: " + "\n" + str(todo) + "\n")


# Translation: Read in the translation from an excel-file and prepare a function to replace german with english
translation = pd.read_excel("translation.xlsx")
translation.columns = ["german", "english"]
def translate(df):
    for deutsch, english in translation.itertuples(index=False,name=None):
        df = df.replace(deutsch,english,regex=True)
    return df

# In[3]: main part: extract the metadata for each excelfile and his sheet(s) in the folder
for excelsheet in todo:

    # read in whole file, then make a list of sheets and select the ones with Rezeptur (can be one or multiple)
    print("\n Working on: "+ os.path.basename(excelsheet))
    excelfile = pd.read_excel(os.path.join(datapath, excelsheet), sheet_name= None) # gives dictionary of dataframes: the keys of the dictionary contain sheet names, and values of the dictionary contain sheet content 
    listofkeys = [i for i in excelfile.keys() if "Rezeptur" in i] # list of Rezeptur dataframes in that dictionary
    print("Relevant sheets in this file: " + str(listofkeys))

    for sheet in listofkeys:

        name = os.path.basename(excelsheet).split(".xl")[0] + " ___ " + sheet
        exceltodf = excelfile[sheet]
        
        # extract header, translate and set proper index & column headers
        exceltodf.iat[17,4] = "Dichte"
        exceltodf = translate(exceltodf)
        specimeninfo = {"Date" : str(exceltodf.columns[9])[:10], 
                "Editor" : exceltodf.iloc[0,10], 
                "Antragsteller" : exceltodf.iloc[2,3],
                "Project number": exceltodf.iloc[3,3],
                "Specimen" : exceltodf.iloc[4,3], 
                "Betonsorte u.-festigkeitsklasse:" : exceltodf.iloc[6,4],
                "Wasserzementwert" : exceltodf.iloc[7,2],
                "Konsistenz" : exceltodf.iloc[7,7],
                "Sieblinie n. DIN 1045" : exceltodf.iloc[8,2],
                "Körnungsziffer" : exceltodf.iloc[8,8],
                "Vol Leim/Zuschlag" : exceltodf.iloc[10,10]
                }
        specimeninfo = pd.Series(specimeninfo)
        exceltodf.rename(columns=exceltodf.iloc[17], inplace = True) # set column header
        exceltodf["Stoffart"] = exceltodf["Stoffart"].str.strip()  # remove WHITESPACE
        
        # manage Zusatzstoff
        howmanyzusatzstoffe = [True if i == "Supplementary cementious materials" else False for i in exceltodf["Stoffart"] ]
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
        df_out.index = df_out.index.map('{0[0]}: {0[1]}'.format)
        df_out = pd.concat([df_out, specimeninfo], ignore_index=False)

         # Write the data to a yaml file
        with open(os.path.join(outputpath, name + '.yaml'), mode='w') as file:
            yaml.dump(df_out.to_dict(), file, sort_keys=False, width=72, indent=4,default_flow_style=None,allow_unicode=True)

            # with open(os.path.join(outputpath, name + '.yaml'), mode='w') as file:
            # yaml.dump(df_out.reset_index().to_dict(orient='records'), file, sort_keys=False, width=72, indent=4,default_flow_style=None,allow_unicode=True)
        
        # Or display it
        # text = yaml.dump(df_out.reset_index().to_dict(orient='records'), sort_keys=False, width=72, indent=4,default_flow_style=None)
        # print(text)
       
print("\n Done.")