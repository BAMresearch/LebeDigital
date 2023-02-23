# Script for the e-module ontology (extracted from CPTO) to map e-module metadata
# by reading every line of that ontology and finding/ # replacing the placeholders.

# import libraries
import yaml
import os
from pathlib import Path
from loguru import logger
import uuid


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
        metadata. Also creates and appends an ID for the specimen.

        Parameter:
        -----
        ontoPath : string
            complete path to Ontology (ttl-format)
        metadataPath : string
            complete path to metadata (yaml-format)
        outputPath : string
            complete path for output

        Output:
        ---
        If no ouput path is given (f.e. for unittesting), the lines will be
        returned. If the "ontoPath" is given for output, the ontology will
        be overwritten. To avoid this, give a new name to create a new ttl-file.

    '''

    # load metadata and get the keys
    metadata = load_metadata(metadataPath)
    keys = list(metadata.keys())  

    # generate ID for the e-module metadata
    specimenID = str(uuid.uuid4())

    # read in the ontology as text linewise, creating a list of lines
    with open(ontoPath, 'r') as file:
        lines = file.readlines()

        # Set up logger
        logger.debug('S T A R T')
        logger.debug('File has ' + str(len(lines)) + ' lines.')
        counter = 0
        usedKeys = [] # to count keys that didn't find a placeholder
        ontoPHcounter = [] # to count all placeholders
        remainingPH = [] # to count the placeholders that recieved no data

        # iterating through the list of lines
        for i in range(len(lines)):

            # create a list of placeholders
            if '_Value$$' in lines[i]:
                ph = lines[i].split("$$")[1]
                ontoPHcounter.append(ph)

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

                # append the specimen-ID name to "key"_    
                key_ = key + "_ "
                if key_ in lines[i]:
                    lines[i] = lines[i].replace(key_, key + "_" + str(specimenID) + " ")

            # append the specimen-ID name to "E-ModulTestSpecimen_"
            if "E-ModulTestSpecimen_ " in lines[i]:
                logger.debug('Found "E-ModulTestSpecimen_" in line ' + str(i + 1) \
                             + ' and appended specimen-ID "' + str(specimenID) + '".')    
                lines[i] = lines[i].replace("E-ModulTestSpecimen_ ", "E-ModulTestSpecimen_" \
                            + str(specimenID) + " ")
                #counter += 1   
                #usedKeys.append(keys[0])

            # depending on shape of specimen, set null pointers            
            if "SpecimenHeight" not in keys:
                if "SpecimenHeight" in lines[i]:
                    placeholder = '$$' + 'SpecimenHeight' + '_Value$$'
                    logger.debug('Specimen shape is cylindrical, set placeholder ' \
                            + placeholder + ' to None.')
                    lines[i] = lines[i].replace(placeholder, str(None))
            if "SpecimenWidth" not in keys:
                if "SpecimenWidth" in lines[i]:
                    placeholder = '$$' + 'SpecimenWidth' + '_Value$$'
                    logger.debug('Specimen shape is cylindrical, set placeholder ' \
                            + placeholder + ' to None.')
                    lines[i] = lines[i].replace(placeholder, str(None))


    ############################ L O G G I N G #############################        

            # create a list of leftover placeholders to see which ones didn't recieve a value
            if '_Value$$' in lines[i]:
                ph = lines[i].split("$$")[1]
                remainingPH.append(ph)

    # for metadata
    unusedKeys = [i for i in keys if i not in usedKeys]
    if len(unusedKeys) > 0:
        logger.warning('Replaced only ' + str(counter) + ' placeholders within the ontology.')
        logger.warning('The following ' + str(len(unusedKeys)) + ' of ' + str(len(keys)) \
                     + ' metadata keys have not been mapped: ')
        logger.warning(unusedKeys)
    else:
        logger.debug('All ' + str(counter) + ' metadata keys have been mapped.')

    # for placeholders
    if len(remainingPH) > 0:
        logger.warning('File has ' + str(len(ontoPHcounter)) + ' placeholders.')
        logger.warning('The following ' + str(len(remainingPH)) + ' of ' + str(len(ontoPHcounter)) \
                    + ' placeholders did not recieve a metadata value: ')
        logger.warning(remainingPH)
    else:
        logger.debug('All ' + str(len(ontoPHcounter)) + ' placeholders within the ontology revieced metadata.')

    ############################ O U T P U T #############################
    if outputPath == None:
        return lines

    else:
        # saving the list again to the file
        with open(outputPath, 'w') as file:
            for line in lines:
                file.write(line)




# # T E M P O R A R Y !!!
# # For my personal testing, will be removed later.

# # defining paths : ONTOLOGY
# ontoDir = Path(__file__).parents[1]
# ontoFile = "../lebedigital/ConcreteOntology/EModuleOntology.ttl"
# ontoPath = os.path.join(ontoDir, ontoFile)

# # defining paths : METADATA
# dataDir = Path(__file__).parents[1]
# dataFile = "../lebedigital/mapping/testMetaData.yaml"  
# dataPath = os.path.join(dataDir, dataFile)

# mappedOntoName = os.path.join(Path(__file__).parents[0], 'EmoduleMappedExmpl.ttl')
# placeholderreplacement(ontoPath, dataPath)#, mappedOntoName)