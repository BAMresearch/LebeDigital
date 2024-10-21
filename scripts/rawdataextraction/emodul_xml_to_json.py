import io
import xml.etree.ElementTree as ET
import json
import datetime
import re


def xml_to_json(xml_binary, mix_json):
    """
    Creates two json-files (E-module and Specimen) containing the
    data from the xml file

    :param xml_binary: xml BLOB from db
    :param mix_json: Json from the Mixture used (binary format)

    :return: [emodul_json, specimen_json]
    """

    xml_file = io.BytesIO(xml_binary)

    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Initialize the JSON object
    emodul_data = {}
    specimen_data = {}

    # Placeholder UUID for emodul_data
    emoduleID = "00000000-0000-4000-8000-000000000000"
    emodul_data['ID'] = emoduleID

    # ID of this specimen, save to specimen metadata and to emodule metadata
    emodul_data['SpecimenID'] = specimen_data['ID'] = emoduleID

    # set experiment lab location to BAM
    emodul_data['Lab'] = 'BAM'

    # Add the raw data path to emodul_data
    emodul_data['RawDataFile'] = "Download"

    # Iterate through each ArrayOfVariableData element
    has_diameter = False
    has_length = False
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
                specimen_data['humanreadableID'] = value
            elif name == 'TestRunName':
                emodul_data['TestRunName'] = value
            elif name == 'E_Modul':
                emodul_data['EModule'] = abs(float(value))
                emodul_data['EModule_Unit'] = var_data.find('Unit').text
            elif name == 'Druckfestigkeit':
                emodul_data['InputCompressiveStrength'] = float(value)
                emodul_data['InputCompressiveStrength_Unit'] = var_data.find('Unit').text
            elif name == 'Durchmesser':
                specimen_data['SpecimenDiameter'] = float(value)
                specimen_data['SpecimenDiameter_Unit'] = var_data.find('Unit').text
                has_diameter = True
            elif name == 'Länge':
                specimen_data['SpecimenLength'] = float(value)
                specimen_data['SpecimenLength_Unit'] = var_data.find('Unit').text
                has_length = True
            elif name == 'Masse':
                specimen_data['SpecimenMass'] = float(value)
                specimen_data['SpecimenMass_Unit'] = var_data.find('Unit').text
            elif name == 'Grundfläche':
                specimen_data['SpecimenBaseArea'] = float(value)
                specimen_data['SpecimenBaseArea_Unit'] = var_data.find('Unit').text
            elif name == 'Rohdichte':
                specimen_data['SpecimenRawDensity'] = float(value)
                specimen_data['SpecimenRawDensity_Unit'] = var_data.find('Unit').text
            elif name == 'Messlänge':
                emodul_data['ExtensometerLength'] = float(value)
                emodul_data['ExtensometerLength_Unit'] = var_data.find('Unit').text
            elif name == 'Dehnung':
                emodul_data['Strain'] = float(value)
                emodul_data['Strain_Unit'] = var_data.find('Unit').text

    # Set specimen shape based on presence of diameter and length
    if has_diameter and has_length:
        specimen_data['SpecimenShape'] = 'Cylinder'
    else:
        specimen_data['SpecimenShape'] = 'Cube'

    json_data = json.loads(mix_json.decode('utf-8'))

    if emodul_data.get("ExperimentDate") and json_data.get("MixingDate"):
        experiment_date = datetime.datetime.strptime(emodul_data['ExperimentDate'],
                                                     '%Y-%m-%dT%H:%M:%S').replace(hour=0, minute=0, second=0)
        mixing_date = datetime.datetime.strptime(json_data["MixingDate"],
                                                 "%Y-%m-%dT%H:%M:%S").replace(hour=0, minute=0, second=0)
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


    return [emodul_data, specimen_data]
