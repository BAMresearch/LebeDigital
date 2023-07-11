import json as js
import os
import numpy as np
import uuid
import datetime
import locale
from loguru import logger
import argparse


def extract_metadata_heat(locationOfRawData, outputDirPath):
    """
        Extracts the metadata from a given datafile (csv).
        Creates also an ID. Returns a dictionary.

        Parameter
        ---------
        locationOfRawData : string
            Path of the csv-file containing the metadata.
        outputDirPath : string
            Path to the directory for output data files (json and csv).

        Output
        -------
        The dict containing the metadata will be returned.
    """

    # create empty dictionary for metadata
    metadata_HoH = {}

    # predefined metadata
    metadata_HoH["Lab"] = "KIT"
    metadata_HoH["ID"] = str(uuid.uuid4())

    metadata_HoH["RawDataFile"] = locationOfRawData
    metadata_HoH["OpenBisRawDataFile"] = None

    # extracted metadata, read in raw data file
    with open(locationOfRawData, "r") as f:
        data = f.readlines()

        # get experiment date and time in protege format YYYY-MM-DDTHH:mm:SS
        locale.setlocale(locale.LC_ALL, 'de_DE.UTF8')
        # starting time
        start = data[2].replace(";", "").replace("\n", "").split(",")[1]
        start = datetime.datetime.strptime(start, "%Y %b %d - %H:%M")
        metadata_HoH["ExperimentStartDate"] = start.strftime('%Y-%m-%dT%H:%M:%S')
        # ending time
        end = data[3].replace(";", "").replace("\n", "").split(",")[1]
        end = datetime.datetime.strptime(end, "%Y %b %d - %H:%M")
        metadata_HoH["ExperimentEndDate"] = end.strftime('%Y-%m-%dT%H:%M:%S')

        mix_names = data[8].replace(";", "").replace("\n", "").split(",")[1:]

        sample_mass = data[9].replace(";", "").replace("\n", "").split(",")[1:]
        sample_mass = [float(m.replace("g", "")) for m in sample_mass]




        # locate numeric experiment data and headers
        # find start line of experiment data
        #print("time;" in data[24].lower()) # check if start is at line 24
        for start_line, line in enumerate(data):
            if "time;" in line.lower():
                headers = line
                break
        #find end line of experiment data
        for end_line, line in enumerate(data):
            if line.startswith(";;") and end_line > start_line:
                break
        headers = headers.split(";")
        numeric_data = []
        logger.debug("Headers: " + str(headers))
        logger.debug("Start line: " + str(start_line))
        logger.debug("End line: " + str(end_line))

        for line in data[start_line+1:end_line]:
            line = line.replace("\n", "").split(";")
            line = [float(l) for l in line]
            numeric_data.append(line)
        numeric_data = np.array(numeric_data)


    #write processed data to json  and csv:
    n = len(mix_names)
    for i, mix in enumerate(mix_names):


        # write metadata to json
        metadata_HoH["mix_name"] = metadata_HoH["humanreadableID"] = mix
        metadata_HoH["sample_mass"] = sample_mass[i]
        mix_headers = headers[:2]
        mix_headers.extend(headers[2+(4*i):(2+4*i)+4])
        metadata_HoH["headers"] = mix_headers
        with open(os.path.join(outputDirPath, mix+"_calo"+".json"), "w") as jsonFile:
            js.dump(metadata_HoH, jsonFile, sort_keys=False, ensure_ascii=False, indent=4)

        # write processed numerical data to csv
        colids = [0, 1]
        colids.extend(range(2+(4*i), (2+4*i)+4))
        np.savetxt(os.path.join(outputDirPath, mix+"_calo"+".csv"), numeric_data[:, colids], delimiter=";", header=";".join(mix_headers).replace(",", " "))

def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script to extract metadata from Heat of Hydration.')
    # input file for raw data
    parser.add_argument('-i', '--input', help='Path to raw data file')
    # output file for metadata json
    parser.add_argument('-o', '--output', help='Path to extracted json files.')
    args = parser.parse_args()

    # default values for testing
    if args.input == None:
        args.input = "../../../usecases/MinimumWorkingExample/Data/HeatOfHydration/calo_example.csv"
    if args.output == None:
        args.output = "../../../usecases/MinimumWorkingExample/Data/HeatOfHydration/processed_data"

    # run extraction and write metadata file
    extract_metadata_heat(args.input, args.output)

if __name__ == "__main__":
    main()



