# Script for checking mix metadata json for duplicates and generate additional
# placeholders for them and add them to the Knowledge graph Template

from loguru import logger
import json
import re
import os

#json_files_directory = '../../usecases/MinimumWorkingExample/mixture/metadata_json_files/'
#json_files = os.listdir(json_files_directory)
json_file_path = '../../usecases/MinimumWorkingExample/mixture/metadata_json_files/2014_08_05 Rezeptur_MI.json'
kg_template_path = '../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl'

def check_duplicate_keys(data, key_prefix):
    # Define the pattern for the given key prefix
    pattern = re.compile(fr'^{re.escape(key_prefix)}\d+_')

    # Find keys matching the pattern
    matching_keys = [key for key in data.keys() if pattern.match(key)]

    # Check if there are multiple keys matching the pattern
    if len(matching_keys) > 1:
        logger.info(f'There are multiple keys following the "{key_prefix}" pattern in the data.')
        logger.info(f'List of matching keys for "{key_prefix}":')
        for key in matching_keys:
            logger.info(f'{key}: {data[key]}')
        return matching_keys
    else:
        logger.info(f'There are no multiple keys following the "{key_prefix}" pattern in the data.')
        return []

def extract_common_prefix(keys):
    common_prefixes = set()

    for key in keys:
        match = re.match(r'^([^_]+_).*', key)
        if match:
            common_prefixes.add(match.group(1))

    return list(common_prefixes)

# Load the JSON data from your file
def check_mix_metadata(json_file_path, kg_template_path):
    #json_files = os.listdir(json_files_directory)
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    #for json_file in json_files:
        #json_file_path = os.path.join(json_files_directory, json_file)
        #with open(json_file_path, 'r') as file:
           # data = json.load(file)

        # Specify the list of key prefixes you're interested in
        key_prefixes_to_check = ['Addition', 'Aggregate', 'Cement', 'Admixture']

        # Initialize the list to store all matching keys
        all_matching_keys = []

        # Call the function for each key prefix
        for key_prefix in key_prefixes_to_check:
            matching_keys = check_duplicate_keys(data, key_prefix)
            all_matching_keys.extend(matching_keys)

        common_prefixes = extract_common_prefix(all_matching_keys)

        # Read the Turtle file line by line
        with open(kg_template_path, 'r') as file:
            turtle_lines = file.readlines()

        # Initialize a set to keep track of existing keys in the TTL file
        existing_keys = set()

        # Iterate over each line and modify subjects
        for i in range(len(turtle_lines)):
            line = turtle_lines[i]
            for prefix in key_prefixes_to_check:
                for common_prefix in common_prefixes:
                    if re.match(f'^{re.escape(prefix)}(\d+)_', common_prefix) and (f'{prefix}_') in line:
                        # Replace the prefix in the line with a value from common_prefixes
                        new_line = line.replace(f'{prefix}_', f'{common_prefix}')
                        turtle_lines[i] = new_line
                        existing_keys.add(common_prefix)
                        break  # Move to the next line after replacement

        # Generate missing blocks for each key in common_prefixes
        for common_prefix in common_prefixes:
            for prefix in key_prefixes_to_check:
                if common_prefix.startswith(f'{prefix}'):
                    # Check if the block is missing for the current common_prefix
                    if common_prefix not in existing_keys:
                        # Generate the block for the missing key
                        new_block = f'\n\nns1:{common_prefix} a owl:NamedIndividual,\n' \
                                    f'        co:BaseMaterial ;\n' \
                                    f'    co:characteristic ns1:{common_prefix}Content_,\n' \
                                    f'        ns1:{common_prefix}Density_ ;\n' \
                                    f'    co:composedOf ns1:{common_prefix}Type_ .\n\n' \
                                    f'ns1:{common_prefix}Content_ a owl:NamedIndividual,\n' \
                                    f'        ns1:Content ;\n' \
                                    f'    co:unit <https://w3id.org/cpto/##{common_prefix}Content_Unit##> ;\n' \
                                    f'    co:value "$${common_prefix}Content_Value$$"^^xsd:float .\n\n' \
                                    f'ns1:{common_prefix}Density_ a owl:NamedIndividual,\n' \
                                    f'        ns1:RelativeDensity ;\n' \
                                    f'    co:unit <https://w3id.org/cpto/##{common_prefix}Density_Unit##> ;\n' \
                                    f'    co:value "$${common_prefix}Density_Value$$"^^xsd:float .\n\n' \
                                    f'ns1:{common_prefix}Type_ a owl:NamedIndividual,\n' \
                                    f'        ns1:{prefix} ;\n' \
                                    f'    co:value "$${common_prefix}Type_Value$$"^^xsd:string .\n\n'

                        # Append the generated block to the modified TTL lines
                        turtle_lines.append(new_block)
                        # Search for the specific line and add missing common prefixes
                    for i in range(len(turtle_lines)):
                        line = turtle_lines[i]
                        if 'co:composedOf' + ' ' + f'ns1:{common_prefix},' in line:
                            # Modify the line to add missing common prefixes
                            new_composed_of = ',\n'.join(
                                f'        ns1:{common_prefix},\n' for common_prefix in common_prefixes if
                                common_prefix not in existing_keys)
                            turtle_lines.insert(i + 1, new_composed_of)
                        # Add the new key to existing_keys
                        existing_keys.add(common_prefix)

        # Join the modified lines back into a string
        modified_turtle_data = ''.join(turtle_lines)


        # Print the modified Turtle data
        #print(modified_turtle_data)

# Open the new TTL file in write mode and write the modified data
    with open(kg_template_path, 'w') as new_file:
        new_file.write(modified_turtle_data)

# Print a message indicating the successful modification and the path to the new TTL file
#logger.info(f'Turtle data has been successfully modified and saved to {turtle_file_path}')


