# Updates:
# This version can accept also excel-files with multiple Rezeptur sheets.

# To-Do: Include xls files. (make sure to have: pip install xlrd) 
#        Avoid any comma in the output.

# ATTENTION: Doesn't work yet with the most recently added mixture datasets from 16.03.2022!

#----------------------------------------------------------------------------------------------------------

# Read excel with pandas, display with iPython, collect files with glob
import pandas as pd
from IPython.display import display
import glob as g
import os

# Locate the mixture data (concrete > Data > Mischungen)
datapath = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', 'Data', 'Mischungen'))

# Displaying setting
pd.set_option('display.max_rows', None)

# Collect a list of excel-Files in the current folder
excelsheetlist = g.glob(datapath + '/*xlsx') # + g.glob('./*xls')  Implement later, causes issues right now
print("There are " + str(len(excelsheetlist)) + " excel-sheets to extract metadata from.")
print(excelsheetlist)

# Translation: Read in the translation from an excel-file and prepare a function to replace german with english
translation = pd.read_excel("translation.xlsx")
translation.columns = ["german", "english"]
def translate(df):
    for deutsch, english in translation.itertuples(index=False,name=None):
        df = df.replace(deutsch,english,regex=True)
    return df
    
# Set up lists to fill inside the loop and build final dataframe from
list_of_dfs = []
list_of_excelsheets = []

# main part: extract the metadata for each excelfile and his sheet(s) in the folder
for excelsheet in excelsheetlist:

    # read in whole file, then make a list of sheets and select the ones with Rezeptur (can be one or multiple)
    print("\n Working on: "+ os.path.basename(excelsheet))
    excelfile = pd.read_excel(excelsheet, sheet_name= None) # gives dictionary of dataframes: the keys of the dictionary contain sheet names, and values of the dictionary contain sheet content 
    listofkeys = [i for i in excelfile.keys() if "Rezeptur" in i] # list of Rezeptur dataframes in that dictionary
    print("Relevant sheets in this file: " + str(listofkeys))

    for sheet in listofkeys:

        list_of_excelsheets.append(os.path.basename(excelsheet) + " : " + sheet)

        # load the sheet and set proper index & column headers & translate
        exceltodf = excelfile[sheet]
        exceltodf.iat[17,2] += " [kg/m^3]"
        exceltodf.iat[17,4] = "Dichte [kg/dm^3]"
        exceltodf = translate(exceltodf)
        exceltodf.rename(columns=exceltodf.iloc[17], inplace = True)
     
        # create new dataframe with only relevant data and chose Index column
        relevant_data = exceltodf.iloc[20:,[0,2,4,8]]
        relevant_data = relevant_data.set_index("Stoffart")
        
        # replace NaN and replace commas for decimal numbers
        relevant_data[relevant_data.columns[-1]] = relevant_data[relevant_data.columns[-1]].fillna("---")
        #relevant_data = .......

        # merge df into a series and add the Zusatzstoff
        df = relevant_data
        df_out = df.stack()
        df_out.index = df_out.index.map('{0[0]}: {0[1]}'.format)
        df_out = pd.concat([df_out, pd.Series([exceltodf.iloc[24,1]],index=["Zusatzstoff"])], ignore_index=False)

        # export information out of loop through list
        list_of_dfs.append(df_out)
       
# concat the dfs from the loop into one
large_df = pd.concat(list_of_dfs, ignore_index=False, axis=1)
large_df.columns = list_of_excelsheets
large_df = large_df.T  # transpose the dataframe (swap columns and rows)

# Display or save as xlsx or csv file
large_df.to_excel("metadata_output.xlsx") 
#large_df.to_csv("metadata_output.csv")
#display(large_df.to_string())





