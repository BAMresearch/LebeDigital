from query import input_emodul_data_for_calibration
import pandas as pd

nameOfExperiment = 'Werner 7.0 M IV E-Modul 28d v. 06.08.14 Probe 6'
df = pd.read_csv(input_emodul_data_for_calibration(nameOfExperiment)['processedDataPath'])

print(df.head())
print('specimen mass: ', input_emodul_data_for_calibration(nameOfExperiment)['specimenMass'])
print('specimen diameter: ', input_emodul_data_for_calibration(nameOfExperiment)['specimenDiameter'])
print('specimen length: ', input_emodul_data_for_calibration(nameOfExperiment)['specimenLength'])