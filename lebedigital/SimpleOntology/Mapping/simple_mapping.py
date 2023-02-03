# Script for the simple ontology to map synthetic data for three specimen to the
# SimpleOntology.ttl file by reading every line of that file and finding/
# replacing the placeholders.

# import libraries
import yaml
import os
from pathlib import Path
from loguru import logger


def load_metadata(dataPath):
    '''
        Load metadata from a given path and return it as dictionary.
        dataPath : string
            Path to the metadata yaml-file.

    '''

    with open(dataPath, 'r') as file:
        try:
            metadata = yaml.safe_load(file)
            return metadata
        except Exception as e:
            logger.error("Path error: " + str(e))


def placeholderreplacement(
        ontoPath,
        metadataPath,
        outputPath=None
):
    '''
        Maps the values of one given metadata file (for one specimen or
        experiment) to a given ontology, by searching within ontology linewise
        for all metadata keys and replacing placeholders with values from the
        metadata. Also appends the name of the specimen.
        Parameter:
        -----
        ontoPath : string
            complete Path to Ontology (ttl-format)
        metadataPath : string
            complete Path to metadata
        outputPath : string
            complete Path for output
        Output:
        ---
        If no ouput path is given (f.e. for unittesting), the lines will be
        returned. If the "ontoPath" is given for output, the ontology will
        be overwritten. To avoid this, give a new name to create a new ttl-file.

    '''

    # load metadata and get the keys
    metadata = load_metadata(metadataPath)
    keys = list(metadata.keys())  # gives ['name', 'ID', 'length', 'diameter', 'shape']

    # read in the ontology as text linewise, creating a list of lines
    with open(ontoPath, 'r') as file:
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
                                 + key + '" with value "' + str(metadata[key]) + '".')
                    lines[i] = lines[i].replace(placeholder, str(metadata[key]))
                    counter += 1
                    usedKeys.append(key)

                # append the specimen name to "key"_
                key_ = key + "_ "
                if key_ in lines[i]:
                    lines[i] = lines[i].replace(key_, key + "_" + str(metadata[keys[0]]) + " ")

            # append the specimen name to "Probe_"
            if "Probe_ " in lines[i]:
                logger.debug('Found "Probe_" in line ' + str(i + 1) \
                             + ' and appended specimen-name "' + str(metadata[keys[0]]) + '".')
                lines[i] = lines[i].replace("Probe_ ", "Probe_" + str(metadata[keys[0]]) + " ")
                counter += 1
                usedKeys.append(keys[0])


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

# defining paths : ONTOLOGY
ontoDir = Path(__file__).parents[1]
ontoFile = "SimpleOntology.ttl"
ontoPath = os.path.join(ontoDir, 'SimpleOntology.ttl')

# defining paths : METADATA
dataDir = Path(__file__).parents[0]
dataFile = "probe1.yaml"  # for now only working with this one, more later
dataPath = os.path.join(dataDir, dataFile)

mappedOntoName = os.path.join(ontoDir, 'SimpleOntology_mappedProbe1.ttl')
placeholderreplacement(ontoPath, dataPath, mappedOntoName)