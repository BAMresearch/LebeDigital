import pandas as pd
import os
import numpy as np
from pathlib import Path


def convert_string_to_number(listStrings):
    listNumbers = []
    for i in listStrings:
        listNumbers.append(float(i.replace(',', '.')))
    return listNumbers


def processed_rawdata(blob):
    """
    Extracts relevant information from the raw data files and returns a dataframe
    Args:
        blob : file from database
        
    Returns:
        DataFrame : processed file
    """
    # Write the blob data to a .dat file
    filename = 'specimen.dat'
    with open(filename, 'wb') as file:
        file.write(blob)

    # Open the .dat file in read mode
    with open(filename, 'r', encoding="utf8", errors='ignore') as file:
        lines = file.readlines()

    #print(len(lines))
    emptyLineIndex = []
    for lineIndex in range(len(lines)):
        if len(lines[lineIndex]) == 1:
            emptyLineIndex.append(lineIndex)

    columns = lines[emptyLineIndex[1] + 2].split('\t')

    rawDataValue = []
    # Calculate the last line index
    last_line_index = len(lines) - 1
    for i in emptyLineIndex:
        if len(emptyLineIndex) == 3:
            for i in range(emptyLineIndex[2] + 4, last_line_index + 1):
                values = convert_string_to_number(lines[i].replace('\n', '').split('\t'))
                if values[3] != 0:
                    rawDataValue.append(values)
    for i in emptyLineIndex:
        if len(emptyLineIndex) > 3:
            if i == emptyLineIndex[1]:
                # Iterate over the first range
                for j in range(i + 4, emptyLineIndex[2]):
                    values = convert_string_to_number(lines[j].replace('\n', '').split('\t'))

                    # Assuming you want to check the fourth column for '0' integers
                    if values[3] != 0:
                        rawDataValue.append(values)
            elif i == emptyLineIndex[3]:
                # Iterate over the second range
                for j in range(i + 4, last_line_index + 1):
                    values = convert_string_to_number(lines[j].replace('\n', '').split('\t'))

                    # Assuming you want to check the fourth column for '0' integers
                    if values[3] != 0:
                        rawDataValue.append(values)



    rawDataDataFrame = pd.DataFrame(columns=[str(n) for n in range(1, 6)], data=rawDataValue)

    processedDataDataFrame = pd.DataFrame(columns=['Force [kN]'
                                                   ],
                                          data=rawDataDataFrame[['5']].values
                                          )
    
    # Clean up the temporary file 
    os.remove(filename)

    return processedDataDataFrame
    #processedDataDataFrame.to_csv(locationOfProcessedData, index=False)




#processed_rawdata('../../../usecases/MinimumWorkingExample/Data/Druckfestigkeit_BAM/20240220_7188_M01/Druckfestigkeiten_7Tage/20240220_7188_M01_W01', '../../../usecases/MinimumWorkingExample/Druckfestigkeit/processeddata')
