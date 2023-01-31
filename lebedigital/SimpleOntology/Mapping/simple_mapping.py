# Script for the simple ontology to map synthetic data for three specimen to the 
# SimpleOntology.ttl file by reading every line of that file and finding/ 
# replacing the placeholders.

# import libraries 
import yaml
import os
from pathlib import Path
from loguru import logger 

###################### LOAD METADATA & ONTOLOGY ###############################

# METADATA
# defining paths
dataDir = Path(__file__).parents[0]
dataFile = "probe1.yaml" #for now only working with this one, more later
dataPath = os.path.join(dataDir, dataFile)

# open yaml-file
with open(dataPath, 'r') as file:
        metadata1 = yaml.safe_load(file)

# get the keys of the metadata
keys = list(metadata1.keys()) # gives ['name', 'ID', 'length', 'diameter', 'shape']
key_append = [(k + "_") for k in keys] # gives ['name_', 'ID_', 'Length_', 'Diameter_', 'shape_']


# ONTOLOGY
# defining paths 
ontoDir = Path(__file__).parents[1]
ontoFile = "SimpleOntology.ttl"
ontoPath = os.path.join(ontoDir,'SimpleOntology.ttl')

# test = "ns2:Diameter_ a ns2:Diameter"
# for i in key_append:
#     if (i+ " ") in test:
#         print(True)

############################ REPLACEMENT-FUNCTION #############################

def placeholderreplacement(
    ontoname,
    outputPath = None
    ):

    '''
        Loads file, searches file linewise for all metadata keys and replaces 
        placeholders with values from the metadata. Also appends the name of the
        specimen to 

        Parameter:
        -----
        ontoname : string
            complete Path to Ontology (ttl-format)
        metadataname : string 
            complete Path to metadata (yaml-format)

        Output:
        ---
        If no ouput path is given (f.e. for unittesting), the lines will be 
        returned. If the "ontoname" is given for output, the ontology will
        be overwritten. To avoid this, give a new name to create a new ttl-file.
    
    '''
    

    # read in the ontology as text linewise, creating a list of lines
    with open(ontoname, 'r') as file:
        lines = file.readlines() 

        # Set up logger
        logger.debug('S T A R T')
        logger.debug('File has ' + str(len(lines)) + ' lines.')
        counter = 0
        usedKeys = []

        # iterating through the list of lines 
        for i in range(len(lines)): 

            # iterate through list of metadata-keys
            for key in keys:

                # create placeholder from keyname f.e. "Length_Value"^^xsd:float
                placeholder = '$$' + key + '_Value$$'

                # if placeholder is in line, replace it with metadata
                if placeholder in lines[i]:
                    logger.debug('Found placeholder "' + placeholder + '" for key "' \
                                 + key + '" with value "' + str(metadata1[key]) + '".')
                    lines[i] = lines[i].replace(placeholder, str(metadata1[key]))
                    counter += 1
                    usedKeys.append(key)

                # append the specimen name to "key"_ 
                key_ = key + "_ "
                if key_ in lines[i]:
                    
                    lines[i] = lines[i].replace(key_, key + "_" + str(metadata1[keys[0]]) + " ")

            # append the specimen name to "Probe_"
            if "Probe_ " in lines[i]:
                logger.debug('Found "Probe_" in line ' + str(i+1) \
                    + ' and appended specimen-name "' + str(metadata1[keys[0]]) + '".')
                lines[i] = lines[i].replace("Probe_ ", "Probe_" + str(metadata1[keys[0]]) + " ")
                counter += 1
                usedKeys.append(keys[0])

            # Special case: shape is replaced by cylindrical and missing underscore
            elif str(metadata1[keys[4]]) in lines[i]:
                lines[i] = lines[i].replace(
                    str(metadata1[keys[4]]), str(metadata1[keys[4]]) + str(metadata1[keys[0]]))


    ############################ L O G G I N G #############################

    unusedKeys = [i for i in keys if i not in usedKeys]
    logger.debug('Replaced ' + str(counter) + ' placeholders within the ontology.')
    if len(unusedKeys) > 0:
        logger.debug('The following ' + str(len(unusedKeys)) + ' of ' + str(len(keys)) \
            + ' metadata keys have not been mapped: ')
        logger.debug(unusedKeys)
    else:
        logger.debug('All metadata keys have been mapped.')

    ############################ O U T P U T #############################
    if outputPath == None:
        return lines

    else:
        # saving the list again to the file
        with open(outputPath, 'w') as file:
            for line in lines:
                file.write(line)


# CREATE EXAMPLE:

mappedOntoName = os.path.join(ontoDir,'SimpleOntology_mappedByScript.ttl')
placeholderreplacement(ontoPath,mappedOntoName)
        



