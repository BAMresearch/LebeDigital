import pandas as pd
import os
from pathlib import Path

def convert_string_to_number(listStrings):
    listNumbers = []
    for i in listStrings:
        listNumbers.append(float(i.replace(',','.')))
    return listNumbers

def processed_data_from_rawdata(locationOfRawData, locationOfProcessedData):
    data = open(os.path.join(locationOfRawData,'specimen.dat'),encoding="utf8", errors='ignore')
    experimentName = os.path.basename(locationOfRawData)
    lines = data.readlines()
                
    emptyLineIndex = []
    for lineIndex in range(len(lines)):
        if len(lines[lineIndex]) == 1:
            emptyLineIndex.append(lineIndex)
    
    columns = lines[emptyLineIndex[1]+2].split('\t')
                
    rawDataValue = []
    for i in range(emptyLineIndex[1] + 4 , emptyLineIndex[2]):
        rawDataValue.append(convert_string_to_number(lines[i].replace('\n','').split('\t')))
                
    rawDataDataFrame = pd.DataFrame(columns=[str(n) for n in range(1,9)], data=rawDataValue)

    processedDataDataFrame = pd.DataFrame(columns= ['Force [kN]',
                                                'Transducer 1[mm]',
                                                'Transducer 2[mm]',
                                                'Transducer 3[mm]'
                                               ],
                                      data = rawDataDataFrame[['4', '6', '7', '8']].values
                                     )
    processedDataDataFrame.to_csv(os.path.join(locationOfProcessedData,'processed_'+ experimentName), index=False)
        
processed_data_from_rawdata('C:\\Users\\vdo\\Desktop\\LeBeDigital\\Code\\minimum_working_example\\ModelCalibration\\usecases\\Concrete\\Example\\Data\\E-modul\\BA Los M V-4', 'C:\\Users\\vdo\\Desktop\\LeBeDigital\\Code\\minimum_working_example\\ModelCalibration\\usecases\\Concrete\\Example\\emodul\\processeddata')