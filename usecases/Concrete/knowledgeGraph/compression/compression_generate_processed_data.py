import pandas as pd
import os
from pathlib import Path


baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
dataFolder = os.path.join(baseDir2,'Data/Druckfestigkeit')
compressionOutputRawData = os.path.join(os.path.join(baseDir0,'compression-processed-data'),'rawdata')
compressionOutputProcessedData = os.path.join(os.path.join(baseDir0,'compression-processed-data'),'processeddata')

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
            
    columns = lines[emptyLineIndex[-1]+2].split('\t')
    numberOfColumn = len(columns)

    rawDataValue = []
    if len(emptyLineIndex) == 3:  
        
        for i in range(emptyLineIndex[2] + 4 , len(lines)):
            rawDataValue.append(convert_string_to_number(lines[i].replace('\n','').split('\t')))
    elif len(emptyLineIndex) == 4:
        for i in range(emptyLineIndex[1] + 4 , emptyLineIndex[2]):
            rawDataValue.append(convert_string_to_number(lines[i].replace('\n','').split('\t')))
        for i in range(emptyLineIndex[3] + 4 , len(lines)):
            rawDataValue.append(convert_string_to_number(lines[i].replace('\n','').split('\t')))
    else:
        print('not a right data format')

    df = pd.DataFrame(columns=[str(n) for n in range(1,numberOfColumn+1)], data=rawDataValue)
    df.to_csv(os.path.join(compressionOutputRawData,csvStoredFileName), index=False)

for folder in os.listdir(dataFolder):
    if folder != 'Maack 8.2 Drucklversuch Probe BK 03 B':
        try:
            dataframe_from_dat(os.path.join(dataFolder,folder,'specimen.dat'),folder.replace(' ','_').replace('.','_')+'.csv')
        except:
            print('something wrong!')
    else:
        print('bad data')

def processed_data_from_rawdata(data):
    df = pd.read_csv(os.path.join(compressionOutputRawData,data))
    processedDataFrame = pd.DataFrame(columns= ['time [s]',
                                                'Force [kN]',
                                                'Transducer [mm]',
                                               ],
                                      data = df[['3', '5', '6']].values
                                     )
    processedDataFrame.to_csv(os.path.join(compressionOutputProcessedData,'processed_'+ data), index=False)
        
for data in os.listdir(compressionOutputRawData):
    processed_data_from_rawdata(data)