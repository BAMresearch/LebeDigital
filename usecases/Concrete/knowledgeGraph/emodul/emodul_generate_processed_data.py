import pandas as pd
import os
from pathlib import Path

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
dataFolder = os.path.join(baseDir2,'Data/E-modul')
emodulOutputRawData = os.path.join(baseDir0,'E-modul-processed-data/rawdata')
emodulOutputProcessedData = os.path.join(baseDir0,'E-modul-processed-data/processeddata')

def convert_string_to_number(listStrings):
    listNumbers = []
    for i in listStrings:
        listNumbers.append(float(i.replace(',','.')))
    return listNumbers

def dataframe_from_dat(dataPath, csvStoredFileName):
    data = open(os.path.join(dataPath),encoding="utf8", errors='ignore')
    lines = data.readlines()
                
    emptyLineIndex = []
    for lineIndex in range(len(lines)):
        if len(lines[lineIndex]) == 1:
            emptyLineIndex.append(lineIndex)
    
    columns = lines[emptyLineIndex[1]+2].split('\t')
                
    rawDataValue = []
    for i in range(emptyLineIndex[1] + 4 , emptyLineIndex[2]):
        rawDataValue.append(convert_string_to_number(lines[i].replace('\n','').split('\t')))
                
    df = pd.DataFrame(columns=[str(n) for n in range(1,9)], data=rawDataValue)
    df.to_csv(os.path.join(emodulOutputRawData,csvStoredFileName), index=False)

for folder in os.listdir(dataFolder):
    try:
        dataframe_from_dat(os.path.join(dataFolder,folder,'specimen.dat'),folder.replace(' ','_').replace('.','_')+'.csv')
    except:
        print('something is wrong!')

def processed_data_from_rawdata(data):
    df = pd.read_csv(os.path.join(emodulOutputRawData,data))
    processedDataFrame = pd.DataFrame(columns= ['Force [kN]',
                                                'Transducer 1[mm]',
                                                'Transducer 2[mm]',
                                                'Transducer 3[mm]'
                                               ],
                                      data = df[['4', '6', '7', '8']].values
                                     )
    processedDataFrame.to_csv(os.path.join(emodulOutputProcessedData,'processed_'+ data), index=False)
        
for data in os.listdir(emodulOutputRawData):
    processed_data_from_rawdata(data)

exampleRawDataframe = pd.read_csv(os.path.join(emodulOutputRawData,'BA-Losert_MI_E-Modul_28d_v__04_08_14_Probe_4.csv'))

exampleProcessedDataframe = pd.read_csv(os.path.join(emodulOutputProcessedData,'processed_BA-Losert_MI_E-Modul_28d_v__04_08_14_Probe_4.csv'))

print('-----------------------')
print('print example raw data csv file')
print(exampleRawDataframe.head())
print('-----------------------')

print('-----------------------')
print('print example processed data csv file')
print(exampleProcessedDataframe.head())
print('-----------------------')