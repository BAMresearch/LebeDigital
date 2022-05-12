import re
import json
import os
import sys
import pandas as pd
import pyaml
from pathlib import Path

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
dataPath = os.path.join(baseDir2,'Data/E-modul')


# the function read each line and return metadata as key and value
def get_metadata_in_one_line(line):
    s = re.sub('\t+', '\t', line)
    s = s.replace('\n','\t')
    result = s.split('\t')[:-1]
    return result

# get file extention and base folder of the file
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


# data = eModul_metadata('/home/dung/Desktop/work/material-digital/data/concrete_research_data/ModelCalibration/usecases/Concrete/Data/E-modul/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4', 'specimen.dat')

# print(data)



def column_data(listName, index, firstKey, secondKey):
    data = []
    for i in listName:
        data.append(i[index][firstKey][secondKey])
    return data

def generate_metadata_yaml_file_from_emodul_BAM(locationOfRawData, locationOfMetaData):

    listDataFolders = os.listdir(locationOfRawData)
    metadata = []
    experimentNameIndex = 0
    dataTypeIndex = 1
    operatorIndex = 2
    dataCollectionIndex = 3

    for folder in listDataFolders:
        path = os.path.join(dataPath,folder)
        
        try:
            dictToYaml = eModul_metadata(path, 'specimen.dat')
            with open(os.path.join(baseDir0,'E-modul-processed-data/metadata_yaml_files/' + folder + '.yaml'), 'w') as yamlFile:
                documents = pyaml.dump(dictToYaml, yamlFile)
            metadata.append(eModul_metadata(path, 'specimen.dat'))
        except:
            print('data file below not in the common structure')
            print(path)


    metadataDict = [{
        'experiment raw name': [i[experimentNameIndex]['experimentName'] for i in metadata],
        'data type': [i[dataTypeIndex]['data type'] for i in metadata],
        'operator time': column_data(metadata, operatorIndex, 'Bediener Information', 'Zeit'),
        'operator timestamp': column_data(metadata, operatorIndex, 'Bediener Information', 'Zeitpunkt'),
        'operator date': column_data(metadata, operatorIndex, 'Bediener Information', 'Datum'),
        'tester': column_data(metadata, operatorIndex, 'Bediener Information', 'Prfer'),
        'sample name': column_data(metadata, operatorIndex, 'Bediener Information', 'Probenbezeichnung'),
        'remark': column_data(metadata, operatorIndex, 'Bediener Information', 'Bemerkungen'),
        'weight': column_data(metadata, operatorIndex, 'Bediener Information', 'Masse'),
        'diameter': column_data(metadata, operatorIndex, 'Bediener Information', 'Durchmesser'),
        'length': column_data(metadata, operatorIndex, 'Bediener Information', 'Lnge'),
        'data collection time': column_data(metadata, dataCollectionIndex, 'Datenerfassung', 'Zeit'),
        'data collection timestamp': column_data(metadata, dataCollectionIndex, 'Datenerfassung', 'Zeitpunkt')
    }]


    dataFrame = pd.DataFrame({
        'experiment raw name': [i[experimentNameIndex]['experimentName'] for i in metadata],
        'data type': [i[dataTypeIndex]['data type'] for i in metadata],
        'operator time': column_data(metadata, operatorIndex, 'Bediener Information', 'Zeit'),
        'operator timestamp': column_data(metadata, operatorIndex, 'Bediener Information', 'Zeitpunkt'),
        'operator date': column_data(metadata, operatorIndex, 'Bediener Information', 'Datum'),
        'tester': column_data(metadata, operatorIndex, 'Bediener Information', 'Prfer'),
        'sample name': column_data(metadata, operatorIndex, 'Bediener Information', 'Probenbezeichnung'),
        'remark': column_data(metadata, operatorIndex, 'Bediener Information', 'Bemerkungen'),
        'weight': column_data(metadata, operatorIndex, 'Bediener Information', 'Masse'),
        'diameter': column_data(metadata, operatorIndex, 'Bediener Information', 'Durchmesser'),
        'length': column_data(metadata, operatorIndex, 'Bediener Information', 'Lnge'),
        'data collection time': column_data(metadata, dataCollectionIndex, 'Datenerfassung', 'Zeit'),
        'data collection timestamp': column_data(metadata, dataCollectionIndex, 'Datenerfassung', 'Zeitpunkt')
    })

    dataFrame.to_csv(os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv'), index = False)

    # dataFrame.to_csv(os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv'), index = False)
    with open(locationOfMetaData, 'w') as yamlFile:
        documents = pyaml.dump(metadataDict, yamlFile)

    return metadataDict

