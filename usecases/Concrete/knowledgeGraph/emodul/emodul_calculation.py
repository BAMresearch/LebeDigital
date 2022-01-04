import pandas as pd
import os
from pathlib import Path
import statistics
import numpy as np
from scipy.signal import argrelextrema
import math

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
metadataFilePath = os.path.join(baseDir0, 'E-modul-processed-data/emodul_metadata.csv')
processedDataPath = os.path.join(baseDir0, 'E-modul-processed-data/processeddata')

df_metadata = pd.read_csv(metadataFilePath)

rawDataFileNames = [df_metadata['experiment raw name'][i] for i in df_metadata.index]

def find_third_element_smaller_than_number(l, m):
    count = 0
    for i, v in enumerate(l):
        if v < m:
            count+=1
            if count == 3:
                return i

emodulValues = []
for file in rawDataFileNames:
    processedDataFileName = 'processed_' + file.replace(' ','_').replace('.','_') + '.csv'
    df = pd.read_csv(os.path.join(processedDataPath, processedDataFileName))
    df['Transducer_average'] = [
                                    statistics.mean([df['Transducer 1[mm]'][i],df['Transducer 2[mm]'][i],df['Transducer 3[mm]'][i]]) for i in df.index
                                ]
    force = np.array(df['Force [kN]'])

    localMaxIndex = argrelextrema(force, np.greater)
    localMaxIndex = list(localMaxIndex[0])
    maxForces = [force[i] for i in localMaxIndex]
    meanMaxForces = statistics.mean(maxForces)
    maxIndex = np.argmax(np.array(maxForces) < meanMaxForces) - 1

    localMinIndex = argrelextrema(force, np.less)
    localMinIndex = list(localMinIndex[0])
    minForces = [force[i] for i in localMinIndex]
    meanMinForces = statistics.mean(minForces)
    minIndex = find_third_element_smaller_than_number(minForces, meanMinForces) 

    specimenDiameter = df_metadata.loc[df_metadata['experiment raw name'] == file]['diameter'].values[0]
    specimenDiameter = float(specimenDiameter.replace(',','.'))

    specimenLength = df_metadata.loc[df_metadata['experiment raw name'] == file]['length'].values[0]
    specimenLength = float(specimenLength.replace(',','.'))
    
    specimenSurfaceArea = math.pi * specimenDiameter**2 / 4000000
    deltaSigma = (maxForces[maxIndex] - minForces[minIndex])*1000/specimenSurfaceArea
    epsilon = list(df['Transducer_average'])
    deltaEpsilon = (epsilon[localMaxIndex[maxIndex]] - epsilon[localMinIndex[minIndex]])/specimenLength

    emodul = deltaSigma/deltaEpsilon
    emodulValues.append(emodul/10**9)

df_metadata['emodul [GPa]'] = emodulValues
df_metadata.to_csv(os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv'), index = False)