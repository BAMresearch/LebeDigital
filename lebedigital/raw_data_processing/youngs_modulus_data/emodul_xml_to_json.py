import os
import xml.etree.ElementTree as ET
import json
import datetime
import re
import uuid
from pathlib import Path
def xml_to_json(xml_file, emodul_json_file, specimen_json_file):
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Initialize the JSON object
    emodul_data = {}
    specimen_data = {}
    # Generate a UUID for emodul_data
    emoduleID = str(uuid.uuid4())
    emodul_data['ID'] = emoduleID

    # ID of this specimen, save to specimen metadata and to emodule metadata
    emodul_data['SpecimenID'] = specimen_data['ID'] = emoduleID

    # set experiment lab location to BAM
    emodul_data['Lab'] = 'BAM'

    # path to XML file
    xml_file_name = os.path.basename(xml_file_path)
    # Extracting the name of the XML file without extension
    xml_file_name_without_extension = os.path.splitext(xml_file_name)[0]
    # Extracting the part before the last underscore
    mix_json_file_name = xml_file_name_without_extension.rsplit('_', 2)[0] + ".json"

    # Adjusting the path to the mix JSON file
    mix_json_file_path = os.path.join("../../../usecases/MinimumWorkingExample/mixture/metadata_json_files/",
                                      mix_json_file_name)

    # save Mixdesign ID to specimen metadata
    try:
        with open(mix_json_file_path, "r", encoding="utf8", errors='ignore') as mixjson:
            mixdesign = json.load(mixjson)
            mixtureID = mixdesign['ID']
            specimen_data['MixtureID'] = mixtureID
            # Extract mixing date
            mixing_date_str = mixdesign['MixingDate']
            mixing_date = datetime.datetime.strptime(mixing_date_str, '%Y-%m-%dT%H:%M:%S')
    except FileNotFoundError:
        raise Exception("No mixdesign json-file found! Can't import the ID and save it to the output!")

    # Get the directory of the XML file
    xml_dir = os.path.dirname(xml_file)

    # Create the path for the raw data file by replacing the extension
    rawdata_file_name = os.path.splitext(os.path.basename(xml_file))[0] + ".xml"
    rawdata_file_path = os.path.join(xml_dir, rawdata_file_name)

    # Add the raw data path to emodul_data
    emodul_data['RawDataFile'] = rawdata_file_path

    # Iterate through each ArrayOfVariableData element
    for array_var_data in root.findall('ArrayOfVariableData'):
        # Iterate through each VariableData element
        for var_data in array_var_data.findall('VariableData'):
            name = var_data.find('Name').text
            value_element = var_data.find('.//Value')
            value = value_element.text

        # Map XML element names to JSON keys
            if name == 'TestRunDate':
                emodul_data['ExperimentDate'] = value
                date_only = datetime.datetime.strptime(value.split(" ")[0], '%d.%m.%Y')
                date_protegeformat = date_only.strftime('%Y-%m-%d') + "T" + value.split(" ")[1]
                emodul_data['ExperimentDate'] = str(date_protegeformat)
            elif name == 'Probenname':
                emodul_data['humanreadableID'] = value
            elif name == 'TestRunName':
                emodul_data['TestRunName'] = value
            elif name == 'E_Modul':
                emodul_data['E_Modul'] = float(value)
                emodul_data['E_Modul_Unit'] = var_data.find('Unit').text
            elif name == 'Druckfestigkeit':
                emodul_data['CompressiveStrength'] = float(value)
                emodul_data['CompressiveStrength_Unit'] = var_data.find('Unit').text
            elif name == 'Durchmesser':
                specimen_data['SpecimenDiameter'] = float(value)
                specimen_data['SpecimenDiameter_Unit'] = var_data.find('Unit').text
            elif name == 'Länge':
                specimen_data['SpecimenLength'] = float(value)
                specimen_data['SpecimenLength_Unit'] = var_data.find('Unit').text
            elif name == 'Masse':
                specimen_data['SpecimenMass'] = float(value)
                specimen_data['SpecimenMass_Unit'] = var_data.find('Unit').text
            elif name == 'Grundfläche':
                specimen_data['SpecimenArea'] = float(value)
                specimen_data['SpecimenArea_Unit'] = var_data.find('Unit').text
            elif name == 'Rohdichte':
                specimen_data['SpecimenRawDensity'] = float(value)
                specimen_data['SpecimenRawDensity_Unit'] = var_data.find('Unit').text
            elif name == 'Messlänge':
                emodul_data['ExtensometerLength'] = float(value)
                emodul_data['ExtensometerLength_Unit'] = var_data.find('Unit').text
            elif name == 'Dehnung':
                emodul_data['Strain'] = float(value)
                emodul_data['Strain_Unit'] = var_data.find('Unit').text

    # Calculate specimen age
    if 'ExperimentDate' in emodul_data and mixing_date:
        # Convert dates to midnight
        experiment_date = datetime.datetime.strptime(emodul_data['ExperimentDate'], '%Y-%m-%dT%H:%M:%S').replace(hour=0, minute=0, second=0)
        mixing_date = mixing_date.replace(hour=0, minute=0, second=0)
        # Calculate age
        specimen_age = (experiment_date - mixing_date).days
        emodul_data['SpecimenAge'] = specimen_age
        emodul_data['SpecimenAge_Unit'] = 'day'
    # Replace \u00b2, \u00b3, etc. with ^2, ^3, etc. in the JSON data
    for key, value in emodul_data.items():
        if isinstance(value, str):
            # Decode JSON string to handle Unicode escape sequences
            decoded_value = json.loads(f'"{value}"')

            # Replace \u00b2, \u00b3, etc. with ^2, ^3, etc.
            decoded_value = re.sub(r'²', '^2', decoded_value)
            decoded_value = re.sub(r'³', '^3', decoded_value)

            emodul_data[key] = decoded_value

    # Replace \u00b2, \u00b3, etc. with ^2, ^3, etc. in the JSON data
    for key, value in specimen_data.items():
        if isinstance(value, str):
            # Decode JSON string to handle Unicode escape sequences
            decoded_value = json.loads(f'"{value}"')

            # Replace \u00b2, \u00b3, etc. with ^2, ^3, etc.
            decoded_value = re.sub(r'²', '^2', decoded_value)
            decoded_value = re.sub(r'³', '^3', decoded_value)

            specimen_data[key] = decoded_value

    # Write emodul_data to JSON file
    with open(emodul_json_file, 'w', encoding='utf-8') as f:
        json.dump(emodul_data, f, indent=4)

    with open(specimen_json_file, 'w', encoding='utf-8') as f:
        json.dump(specimen_data, f, indent=4)

#def main():
    # Folder containing XML files
#    xml_folder = "../../../usecases/MinimumWorkingExample/Data/E-Modul_28_Tage"

    # Folder for JSON output files
#    json_folder = "path/to/your/json_folder"

    # Iterate through XML files in the folder
#    for file_name in os.listdir(xml_folder):
#        if file_name.endswith(".xml"):
#            xml_file_path = os.path.join(xml_folder, file_name)
#            json_file_path = os.path.join(json_folder, file_name.replace(".xml", ".json"))
#            xml_to_json(xml_file_path, json_file_path)
#            print("Conversion successful. JSON file saved at:", json_file_path)

#if __name__ == "__main__":
#    main()


# Example usage:
xml_file_path = "../../../usecases/MinimumWorkingExample/Data/E-Modul_28_Tage/20240220_7188_M01_Z04_E.xml"
emodul_json_file_path = "../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/20240220_7188_M01_Z04_E.json"
specimen_json_file_path = "../../../usecases/MinimumWorkingExample/emodul/metadata_json_files/20240220_7188_M01_Z04_E_Specimen.json"

xml_to_json(xml_file_path, emodul_json_file_path, specimen_json_file_path)

print("Conversion successful. JSON files saved at:", emodul_json_file_path, specimen_json_file_path)


