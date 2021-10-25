import re
import json
import os
import sys
import pandas as pd
from pathlib import Path

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
dataPath = os.path.join(baseDir2,'Data/Druckfestigkeit')

csvPath = os.path.join(os.path.join(baseDir0,'compression-processed-data'),'compression_metadata.csv')

listDataFolders = os.listdir(dataPath)

def get_metadata_in_one_line(line):
    s = re.sub('\t+', '\t', line)
    s = s.replace('\n','\t')
    result = s.split('\t')[:-1]
    return result

def get_file_information(path, file):
    completedFileName, fileExtention = os.path.splitext(os.path.join(path,file))
    folderName = os.path.basename(path)
    return fileExtention, folderName

def eModul_metadata(filePath, fileName ):
    
    if sys.platform == 'windows':
        data = open(r"{}\{}".format(filePath,fileName),encoding="utf8", errors='ignore')
    else:
        data = open(r"{}/{}".format(filePath,fileName),encoding="utf8", errors='ignore')
    # get type of the file and name of base folder
    fileExtention, experimentName = get_file_information(filePath, fileName)

    if fileExtention == '.dat':
        fileExtention = 'DATFile'
    elif fileExtention == '.csv':
        fileExtention = 'CSVFile'
    elif fileExtention == '.cad':
        fileExtention = 'CADFile'
    else:
        fileExtention = 'DocumentFile'

    dataType = {
        'data type': fileExtention
    }

    # read data file line by line
    lines = data.readlines()
    # get empty lines (where start and end the header)
    emptyLineIndex = []
    for lineIndex in range(len(lines)):
        if len(lines[lineIndex]) == 1:
            emptyLineIndex.append(lineIndex)

    # service information of the experiment, it should be in between the first two empty lines
    serviceInformation = []
    for ind in range(emptyLineIndex[0]+1,emptyLineIndex[1]):
        serviceInformation.append(get_metadata_in_one_line(lines[ind]))

    # generate service information dictionary
    serviceInformationDict = {'Bediener Information':{
        serviceInformation[0][1].replace(':','').strip():serviceInformation[0][2],
        'Zeitpunkt':serviceInformation[0][4],
    }}
    for i in range(len(serviceInformation)-2):
        
        serviceInformationDict['Bediener Information'][serviceInformation[i+1][0].replace(':','').strip()] = serviceInformation[i+1][1]

    # data collection dictionary
    dataCollectionInformationDict = {
        'Datenerfassung': {
            get_metadata_in_one_line(lines[emptyLineIndex[1]+1])[1].replace(':',''):get_metadata_in_one_line(lines[emptyLineIndex[1]+1])[2],
            'Zeitpunkt':get_metadata_in_one_line(lines[emptyLineIndex[1]+1])[4]
        }
    }

    # columns name and unit

    columnsDict = { 'Column Data': {
        'Columns Name':get_metadata_in_one_line(lines[emptyLineIndex[1]+2]),
        'Column Unit':get_metadata_in_one_line(lines[emptyLineIndex[1]+3])
        }
        
    }
    experimentNameDict = {'experimentName': experimentName}
    # aggregate the metadata
    metadata = [experimentNameDict,dataType, serviceInformationDict,dataCollectionInformationDict,columnsDict]
    # jsonData = json.dumps(serviceInformationDict,ensure_ascii=False)
    
    return metadata

metadata = []
for folder in listDataFolders:
    if folder != 'Maack 8.2 Drucklversuch Probe BK 03 B':
        path = os.path.join(dataPath,folder)
        try:
            metadata.append(eModul_metadata(path, 'specimen.dat'))
        except:
            print('something wrong')
            print(folder)
    else:
        print('not a right data')

def column_data(listName, index, firstKey, secondKey):
    data = []
    for i in listName:
        data.append(i[index][firstKey][secondKey])
    return data

experimentNameIndex = 0
dataTypeIndex = 1
operatorIndex = 2
dataCollectionIndex = 3

dataFrame = pd.DataFrame({
    'experiment raw name': [i[experimentNameIndex]['experimentName'] for i in metadata],
    'data type': [i[dataTypeIndex]['data type'] for i in metadata],
    'tester': column_data(metadata, operatorIndex, 'Bediener Information', 'Prfer'),
    'sample name': [i[experimentNameIndex]['experimentName'] for i in metadata],
    'weight': column_data(metadata, operatorIndex, 'Bediener Information', 'Masse [g]'),
    'diameter': column_data(metadata, operatorIndex, 'Bediener Information', 'Durchmesser [mm]'),
    'height': column_data(metadata, operatorIndex, 'Bediener Information', 'Hhe [mm]'),
    'data collection timestamp': column_data(metadata, dataCollectionIndex, 'Datenerfassung', 'Zeitpunkt')
})

dataFrame.to_csv(csvPath, index=False)