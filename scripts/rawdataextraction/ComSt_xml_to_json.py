import xml.etree.ElementTree as ET
import json
import datetime
import re
import uuid
import io

def comSt_xml_to_json(xml_binary, mix_json):
    """
    Creates two json-files (E-module and Specimen) containing the
    data from the xml file

    :param xml_binary: xml BLOB from db
    :param mix_json: Json from the Mixture used (binary format)

    :return: [ComSt_json, specimen_json]
    """

    xml_file = io.BytesIO(xml_binary)

    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Initialize the JSON object
    ComSt_data = {}
    specimen_data = {}
    # Generate a UUID for emodul_data
    ComStID = str(uuid.uuid4())
    ComSt_data['ID'] = ComStID

    # ID of this specimen, save to specimen metadata and to emodule metadata
    ComSt_data['SpecimenID'] = specimen_data['ID'] = ComStID

    # set experiment lab location to BAM
    ComSt_data['Lab'] = 'BAM'

    # Process mix_json data
    try:
        mixdesign = json.loads(mix_json.decode('utf-8'))
        # Extract mixing date
        mixing_date_str = mixdesign.get("MixingDate")
        mixing_date = datetime.datetime.strptime(mixing_date_str, '%Y-%m-%dT%H:%M:%S')
    except Exception as e:
        raise Exception(f"Error processing mix_json: {str(e)}")

    # Add the raw data path to emodul_data
    ComSt_data['RawDataFile'] = "Download"

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
                date_only = datetime.datetime.strptime(value.split(" ")[0], '%d.%m.%Y')
                date_protegeformat = date_only.strftime('%Y-%m-%d') + "T" + value.split(" ")[1]
                ComSt_data['ExperimentDate'] = str(date_protegeformat)
            elif name == 'Probenname':
                specimen_data['humanreadableID'] = value
            elif name == 'TestRunName':
                ComSt_data['TestRunName'] = value
            elif name == 'Druckfestigkeit':
                ComSt_data['CompressiveStrength'] = float(value)
                ComSt_data['CompressiveStrength_Unit'] = var_data.find('Unit').text
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

    # Set specimen shape based on presence of diameter and length
    if has_diameter and has_length:
        specimen_data['SpecimenShape'] = 'Cylinder'
    else:
        specimen_data['SpecimenShape'] = 'Cube'

    # Calculate specimen age
    if 'ExperimentDate' in ComSt_data and mixing_date:
        # Convert dates to midnight
        experiment_date = datetime.datetime.strptime(ComSt_data['ExperimentDate'], '%Y-%m-%dT%H:%M:%S').replace(hour=0, minute=0, second=0)
        mixing_date = mixing_date.replace(hour=0, minute=0, second=0)
        # Calculate age
        specimen_age = (experiment_date - mixing_date).days
        ComSt_data['SpecimenAge'] = specimen_age
        ComSt_data['SpecimenAge_Unit'] = 'day'

    # Replace special characters in both objects
    for data in [ComSt_data, specimen_data]:
        for key, value in data.items():
            if isinstance(value, str):
                # Decode JSON string to handle Unicode escape sequences
                decoded_value = json.loads(f'"{value}"')
                # Replace \u00b2, \u00b3, etc. with ^2, ^3, etc.
                decoded_value = re.sub(r'²', '^2', decoded_value)
                decoded_value = re.sub(r'³', '^3', decoded_value)
                data[key] = decoded_value

    return [ComSt_data, specimen_data]


