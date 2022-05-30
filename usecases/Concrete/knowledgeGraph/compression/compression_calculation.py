import pandas as pd
import statistics
import numpy as np
from scipy.signal import argrelextrema
import pandas as pd
import os
from pathlib import Path
import math

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
metadataFilePath = os.path.join(baseDir0, 'compression-processed-data/compression_metadata.csv')
processedDataPath = os.path.join(baseDir0, 'compression-processed-data/processeddata')

df_metadata = pd.read_csv(metadataFilePath)
processedDataFileNames = [df_metadata['experiment raw name'][i] for i in df_metadata.index]

compressionStrengthValues = []
for file in processedDataFileNames:
    processedDataFileName = 'processed_' + file.replace(' ','_').replace('.','_') + '.csv'
    df = pd.read_csv(os.path.join(processedDataPath, processedDataFileName))
    forces = [i for i in df['Force [kN]'] if i < 0]
    minIndex = np.where(forces == np.min(forces))
    maxForce = abs(np.min(forces)*1000)

    specimenDiameter = df_metadata.loc[df_metadata['experiment raw name'] == file]['diameter'].values[0]
    specimenDiameter = float(specimenDiameter.replace(',','.'))
    
    specimenSurfaceArea = math.pi * specimenDiameter**2 / 4000000
    compressionStrength = maxForce/specimenSurfaceArea
    compressionStrengthValues.append(compressionStrength/10**6)

df_metadata['compressive_strength [MPa]'] = compressionStrengthValues
df_metadata.to_csv(os.path.join(baseDir0,'compression-processed-data/compression_metadata.csv'), index = False)