# Script for a knowledge graph template derived from an ontology
# to map metadata by reading every line of the template and finding/ replacing the placeholders.
# Logging through loguru, you can ignore "debug" messages. "Warning" appears if not
# everything has been mapped.

# import libraries
import json
import os
from pathlib import Path
from loguru import logger
import uuid
import argparse
from lebedigital.mapping.unit_conversion import unit_conversion
from lebedigital.mapping.check_duplicate import check_mix_metadata

def load_metadata(dataPath):
    '''
        Load metadata from a given path and return it as dictionary.
        dataPath : string
            Path to the metadata json-file.

    '''

    with open(dataPath, 'r') as file:
        try:
            metadata = json.load(file)
            return metadata
        except Exception as e:
            logger.error("Path error: " + str(e))


def generate_placeholder(key, type="Value"):
    '''
        Generates a placeholder (str), standard (type="Value") in the format $$key_Value$$ for a given key.
        If type="Unit", then the placeholder has the format ##key_Unit##.
        This function should allow to easily change the structure of the placeholder
        given in the template without having to rewrite the function placeholderreplacement.
        Just change the structure here.
    '''

    if type == "Value":
        placeholder = '$$' + str(key) + '_Value$$'
    else:
        placeholder = '##' + str(key) + '##'

    return placeholder


def placeholderreplacement(kgPath, metadataPath):
    '''
        Maps the values of one given metadata file (for one specimen or
        experiment) to a given knowledge graph (KG) template, by searching through
        it linewise for all metadata keys and replacing placeholders with values 
        from the metadata. Also appends an ID for the specimen.

        Parameter:
        -----
        kgPath : string
            complete path to KG template with placeholders (ttl-format)
        metadataPath : string
            complete path to metadata (json-format)

        Output:
        ---
        Returns list of lines that consist of the template with mapped metadata.

    '''

    # load metadata, convert the units through module and get the keys
    metadata = load_metadata(metadataPath)
    metadata = unit_conversion(metadata)
    keys = list(metadata.keys())
    keys_unit = [i for i in keys if "_Unit" in i]
    keys_value = [i for i in keys if i not in keys_unit]

    # search metadata for duplicates ...

    # duplicate template part regarding the key

    # import ID from metadata to append to all instances
    metadataID = metadata["ID"]

    # read in the KG template as text linewise, creating a list of lines
    with open(kgPath, 'r') as file:
        lines = file.readlines()

        # Set up logger
        logger.debug('S T A R T')
        logger.debug('Loaded ttl-File has ' + str(len(lines)) + ' lines.')
        usedKeys = []  # to count keys that found a placeholder
        kgPHcounter = []  # to count all placeholders
        remainingPH = []  # to count the placeholders that recieved no data

        # iterating through the list of lines
        for i in range(len(lines)):

            # create a list of placeholders
            if '_Value$$' in lines[i]:
                ph = lines[i].split("$$")[1]
                kgPHcounter.append(ph)
            if '_Unit##' in lines[i]:
                ph = lines[i].split("##")[1]
                kgPHcounter.append(ph)

            # iterate through list of metadata-keys
            for key in keys_value:

                placeholder = generate_placeholder(key)

                # if placeholder is in line, replace it with metadata
                if placeholder in lines[i]:
                    logger.debug('Found value placeholder "' + placeholder + '" for key "' \
                                 + key + '" with value "' + str(metadata[key]) + '".')
                    lines[i] = lines[i].replace(placeholder, str(metadata[key]))
                    usedKeys.append(key)

                # append the specimen-ID name to "key"_ , works for most keys, except
                # some keys below
                key_ = key + "_ "
                if key_ in lines[i]:
                    lines[i] = lines[i].replace(key_, key + "_" + str(metadataID) + " ")

            # iterate through list of unit-keys
            for key in keys_unit:

                # if unit is in line, replace unit-placeholder with proper unit
                placeholder_unit = generate_placeholder(key, "Unit")

                # if placeholder is in line, replace it with unit
                if placeholder_unit in lines[i]:
                    logger.debug('Found unit placeholder "' + placeholder_unit + '" for key "' \
                                 + key + '" with value "' + str(metadata[key]) + '".')
                    replace = lines[i].replace(">", "<")
                    replace = replace.split("<")[1]
                    lines[i] = lines[i].replace(replace, str(metadata[key]))
                    usedKeys.append(key)


            # append the specimen-ID name to the exceptions 
            if "_," in lines[i]:
                #logger.debug('Appended specimen-ID in line ' + str(i + 1) \
                #             + ' to ' + str(lines[i].split("_,")[0] + "_,") + '".')    
                lines[i] = lines[i].replace("_,", "_" + str(metadataID) + ",")
            if "_ " in lines[i]:
                #logger.debug('Appended specimen-ID in line ' + str(i + 1) \
                #             + ' to ' + str(lines[i].split("_ ")[0] + "_ ") + '".')    
                lines[i] = lines[i].replace("_ ", "_" + str(metadataID) + " ")


    ############################ L O G G I N G #############################        

            # create a list of leftover placeholders to see which ones didn't receive a value
            if '_Value$$' in lines[i]:
                ph = lines[i].split("$$")[1]
                remainingPH.append(ph)
                lines[i] = ''
                lines[i - 1] = lines[i - 1].replace(';', '.')
            elif '_Unit##' in lines[i]:
                ph = lines[i].split("##")[1]
                remainingPH.append(ph)
                lines[i] = ''
                lines[i + 1] = ''
                lines[i - 1] = lines[i - 1].replace(';', '.')
            elif '"None"^^xsd' in lines[i]:
                lines[i] = ''
                lines[i - 1] = lines[i - 1].replace(';', '.')
                #lines[i - 1] = ''
                #lines[i - 2] = ''
            elif '<None>' in lines[i]:
                lines[i] = ''
                lines[i - 1] = ''
                lines[i - 2] = ''

    # for metadata
    unusedKeys = [i for i in keys if i not in usedKeys]
    if len(unusedKeys) > 0:
        logger.warning('Mapped ' + str(len(usedKeys)) + ' keys to the KG template.')
        logger.warning(usedKeys)
        logger.warning('The following ' + str(len(unusedKeys)) + ' of ' + str(len(keys)) \
                     + ' metadata keys have not been mapped: ')
        logger.warning(unusedKeys)
    else:
        logger.debug('All ' + str(len(usedKeys)) + ' metadata keys have been mapped.')

    # for placeholders
    if len(remainingPH) > 0:
        logger.warning('File has ' + str(len(kgPHcounter)) + ' placeholders.')
        logger.warning('The following ' + str(len(list(set(remainingPH)))) + ' of ' + str(len(kgPHcounter)) \
                    + ' placeholders did not receive a metadata value: ')
        logger.warning(list(set(remainingPH)))
    else:
        logger.debug('All ' + str(len(kgPHcounter)) + ' placeholders within the KG template received metadata.')

    return lines


def mapping(KGtemplatePath, metadataPath, outputPath):
    """Returns a ttl-file based on a KG template, where placeholders have been replaced with data from a given metadata file

    Parameters
    ----------
    metadataPath : string
        Path to the metadata file
    KGtemplatePath : string
        Path to the KG template (ttl-file)
    outputPath : string
        Path to where the mapped KG template should be stored
    """

    # mapping
    mappedKG = placeholderreplacement(KGtemplatePath, metadataPath)

    # writing to ttl file
    with open(outputPath, 'w') as file:
        for line in mappedKG:
            file.write(line)

def main():
    # create parser
    parser = argparse.ArgumentParser(description='Script for mapping metadata to a Knowledge graph template.')
    # input file for KGs and metadata
    parser.add_argument('-i', '--input', help='Paths to knowledge graph template and to metadata')
    # output for mapped KGs
    parser.add_argument('-o', '--output', help='Path to the mapped graph.')
    args = parser.parse_args()

    # default values for testing of my script
    if args.input == None:
        args.input = ['../../lebedigital/ConcreteOntology/CompressiveStrength_KG_Template.ttl',
                     '../../usecases/demonstrator/KIT_Data/CompressiveStrength_json_files/mix0_CompressiveStrength/KIT_21-1605_M1_1d_W1.json',
                      '../../lebedigital/ConcreteOntology/Specimen_KG_Template.ttl',
                      '../../usecases/demonstrator/KIT_Data/CompressiveStrength_json_files/mix0_CompressiveStrength/KIT_CS_M1_1d_W1_Specimen.json',
                      '../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template.ttl',
                      '../../usecases/MinimumWorkingExample/mixture/metadata_json_files/2019_06_26 Klimek Geschossdecke_Quarzkies.json',
                      '../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template_modified.ttl'
                      ]
    if args.output == None:
        args.output = ['../../usecases/MinimumWorkingExample/Mapping_Example/ComStMapped.ttl',
                       '../../usecases/MinimumWorkingExample/Mapping_Example/SpecimenMapped.ttl',
                       '../../usecases/MinimumWorkingExample/Mapping_Example/MixMapped.ttl',
                       '../../lebedigital/ConcreteOntology/MixtureDesign_KG_Template_modified.ttl']

    # Check mix metadata and generate additional placeholders
    check_mix_metadata(args.input[5], args.input[4], args.output[3])

    # run extraction and write metadata file
    #emodule
    #mapping(args.input[0], args.input[1], args.output[0])
    #specimen
    #mapping(args.input[2], args.input[3], args.output[1])
    #mix
    mapping(args.input[6], args.input[5], args.output[2])

if __name__ == "__main__":
    main()
