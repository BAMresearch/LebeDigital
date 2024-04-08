import pandas as pd
import os
import numpy as np
from pathlib import Path


def convert_string_to_number(listStrings):
    listNumbers = []
    for i in listStrings:
        listNumbers.append(float(i.replace(',', '.')))
    return listNumbers


def processed_rawdata(locationOfRawData, locationOfProcessedData):
    """
    Extracts relevant information from the raw data files and stores as a csv file
    Args:
        locationOfRawData (): string
            The path to the raw data
        locationOfProcessedData (): string
            The path where the processed data needs to be stored

    Returns:

    """
    data = open(os.path.join(locationOfRawData, 'specimen.dat'), encoding="utf8", errors='ignore')


    experimentName = os.path.basename(locationOfRawData)
    lines = data.readlines()
    print(len(lines))
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



    rawDataDataFrame = pd.DataFrame(columns=[str(n) for n in range(1, 7)], data=rawDataValue)

    processedDataDataFrame = pd.DataFrame(columns=['Force [kN]'
                                                   ],
                                          data=rawDataDataFrame[['5']].values
                                          )
    processedDataDataFrame.to_csv(locationOfProcessedData, index=False)




processed_rawdata('C:/develop/lebedigital-new/Lebedigital/usecases/MinimumWorkingExample/Data/Druckfestigkeit/Wolf 8.2 Probe 2', 'C:/develop/lebedigital-new/Lebedigital/usecases/MinimumWorkingExample/Druckfestigkeit/processeddata')



