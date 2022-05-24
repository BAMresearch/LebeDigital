import re
import json
import os
import sys
import pandas as pd
import pyaml
from pathlib import Path


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

def eModul_metadata(locationOfRawData, fileName,locationOfMetaData ):
    
    if sys.platform == 'win32':
        data = open(r"{}\{}".format(locationOfRawData,fileName),encoding="utf8", errors='ignore')
    else:
        data = open(r"{}/{}".format(locationOfRawData,fileName),encoding="utf8", errors='ignore')
    # get type of the file and name of base folder
    fileExtention, experimentName = get_file_information(locationOfRawData, fileName)

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

    metadataDict = {
        'experimentName': metadata[0]['experimentName'],
        'dataType': metadata[1]['data type'],
        'operatorTime': metadata[2]['Bediener Information']['Zeit'],
        'operatorTimestamp': metadata[2]['Bediener Information']['Zeitpunkt'],
        'operatorDate': metadata[2]['Bediener Information']['Datum'],
        'tester': metadata[2]['Bediener Information']['Prfer'],
        'specimenName': metadata[2]['Bediener Information']['Probenbezeichnung'],
        'remark': metadata[2]['Bediener Information']['Bemerkungen'],
        'weight': metadata[2]['Bediener Information']['Masse'],
        'diameter': metadata[2]['Bediener Information']['Durchmesser'],
        'length': metadata[2]['Bediener Information']['Lnge'],
    }

    with open(locationOfMetaData, 'w') as yamlFile:
        documents = pyaml.dump(metadataDict, yamlFile)
    return metadataDict



# eModul_metadata('C:\\Users\\vdo\\Desktop\\LeBeDigital\\Code\\minimum_working_example\\ModelCalibration\\usecases\\Concrete\\Example\\Data\\E-modul\\BA Los M V-4', 'specimen.dat', 'C:\\Users\\vdo\\Desktop\\LeBeDigital\\Code\\minimum_working_example\\ModelCalibration\\usecases\\Concrete\\Example\\emodul\\metadata_yaml_files')