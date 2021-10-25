from compression_query import input_compression_data_for_calibration
import pandas as pd

nameOfExperiment = 'Hüsken Probe 1-2'
df = pd.read_csv(input_compression_data_for_calibration(nameOfExperiment)['processedDataPath'])

print(df.head())
print('specimen mass: ', input_compression_data_for_calibration(nameOfExperiment)['specimenMass'])
print('specimen diameter: ', input_compression_data_for_calibration(nameOfExperiment)['specimenDiameter'])
print('specimen height: ', input_compression_data_for_calibration(nameOfExperiment)['specimenHeight'])