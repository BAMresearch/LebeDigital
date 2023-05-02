# Script for a knowledge graph template derived from an ontology
# to map metadata by reading every line of the template and finding/ replacing the placeholders.
# Logging through loguru, you can ignore "debug" messages. "Warning" appears if not
# everything has been mapped.

# import libraries
import yaml
import os
from pathlib import Path
from loguru import logger
import uuid
#from lebedigital.mapping.unit_conversion import unit_conversion


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


def generate_placeholder(key):
    '''
        Generates a placeholder (str) in the format $$key_Value$$ for a given key.
        This function should allow to easily change the structure of the placeholder
        given in the template without having to rewrite the function placeholderreplacement.
        Just change the structure here.
    '''

    placeholder = '$$' + str(key) + '_Value$$'
    return placeholder


def placeholderreplacement(
        kgPath,
        metadataPath,
        outputPath=None
        ):
    '''
        Maps the values of one given metadata file (for one specimen or
        experiment) to a given knowledge graph (KG) template, by searching through
        it linewise for all metadata keys and replacing placeholders with values 
        from the metadata. Also creates and appends an ID for the specimen.

        Parameter:
        -----
        kgPath : string
            complete path to KG template with placeholders (ttl-format)
        metadataPath : string
            complete path to metadata (yaml-format)
        outputPath : string
            complete path for output

        Output:
        ---
        If no output path is given (f.e. for unittesting), the lines will be
        returned. If the "kgPath" is given for output, the KG template will
        be overwritten. To avoid this, give a new name to create a new ttl-file.

    '''

    # load metadata and get the keys
    metadata = load_metadata(metadataPath)
    keys = list(metadata.keys())  
    #print(unit_conversion(metadata))
    # generate ID
    specimenID = str(uuid.uuid4())

    # read in the KG template as text linewise, creating a list of lines
    with open(kgPath, 'r') as file:
        lines = file.readlines()

        # Set up logger
        logger.debug('S T A R T')
        logger.debug('Loaded ttl-File has ' + str(len(lines)) + ' lines.')
        usedKeys = [] # to count keys that found a placeholder
        kgPHcounter = [] # to count all placeholders
        remainingPH = [] # to count the placeholders that recieved no data

        # iterating through the list of lines
        for i in range(len(lines)):

            # create a list of placeholders
            if '_Value$$' in lines[i]:
                ph = lines[i].split("$$")[1]
                kgPHcounter.append(ph)

            # iterate through list of metadata-keys
            for key in keys:

                placeholder = generate_placeholder(key)

                # if placeholder is in line, replace it with metadata
                if placeholder in lines[i]:
                    logger.debug('Found placeholder "' + placeholder + '" for key "' \
                                 + key + '" with value "' + str(metadata[key]) + '".')
                    lines[i] = lines[i].replace(placeholder, str(metadata[key]))
                    usedKeys.append(key)

                # append the specimen-ID name to "key"_ , works for most keys, except 
                # some keys below
                key_ = key + "_ "
                if key_ in lines[i]:
                    lines[i] = lines[i].replace(key_, key + "_" + str(specimenID) + " ")


            # append the specimen-ID name to the exceptions 
            if "_," in lines[i]:
                #logger.debug('Appended specimen-ID in line ' + str(i + 1) \
                #             + ' to ' + str(lines[i].split("_,")[0] + "_,") + '".')    
                lines[i] = lines[i].replace("_,", "_" + str(specimenID) + ",")
            if "_ " in lines[i]:
                #logger.debug('Appended specimen-ID in line ' + str(i + 1) \
                #             + ' to ' + str(lines[i].split("_ ")[0] + "_ ") + '".')    
                lines[i] = lines[i].replace("_ ", "_" + str(specimenID) + " ")


            # ID-key is not given by metadata but created in this script, so map it now:
            if generate_placeholder("SpecimenID") in lines[i]:
                    logger.debug('Found placeholder "' + generate_placeholder("SpecimenID")+ '".')
                    lines[i] = lines[i].replace(generate_placeholder("SpecimenID"), str(specimenID))


    ############################ L O G G I N G #############################        

            # create a list of leftover placeholders to see which ones didn't receive a value
            if '_Value$$' in lines[i]:
                ph = lines[i].split("$$")[1]
                remainingPH.append(ph)

    # for metadata
    unusedKeys = [i for i in keys if i not in usedKeys]
    if len(unusedKeys) > 0:
        logger.warning('Mapped only ' + str(len(usedKeys)) + ' keys to the KG template.')
        logger.warning('The following ' + str(len(unusedKeys)) + ' of ' + str(len(keys)) \
                     + ' metadata keys have not been mapped: ')
        logger.warning(unusedKeys)
    else:
        logger.debug('All ' + str(len(usedKeys)) + ' metadata keys have been mapped.')

    # for placeholders
    if len(remainingPH) > 0:
        logger.warning('File has ' + str(len(kgPHcounter)) + ' placeholders.')
        logger.warning('The following ' + str(len(remainingPH)) + ' of ' + str(len(kgPHcounter)) \
                    + ' placeholders did not receive a metadata value: ')
        logger.warning(remainingPH)
    else:
        logger.debug('All ' + str(len(kgPHcounter)) + ' placeholders within the KG template received metadata.')

    ############################ O U T P U T #############################
    if outputPath == None:
        return lines

    else:
        # saving the list again to the file
        with open(outputPath, 'w') as file:
            for line in lines:
                file.write(line)
