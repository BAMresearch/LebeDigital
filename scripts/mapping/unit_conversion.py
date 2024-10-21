from loguru import logger


def unit_conversion(input_metadata):
    """
    Script to convert the units extracted from the raw data into valid ontology
    instances.

    :param input_metadata: dictionary of metadata
    :return: dictionary of converted units
    """

    # Define the unit mappings as a dictionary
    unit_mappings = {"mm": "http://qudt.org/vocab/unit/MilliM",
                     "g": "http://qudt.org/vocab/unit/GM",
                     "kN": "https://qudt.org/vocab/unit/KN",
                     "kg/m^3": "http://qudt.org/vocab/unit/KiloGM-PER-M3",
                     "kg/dm^3": "http://qudt.org/vocab/unit/KiloGM-PER-DeciM3",
                     "dm^3": "http://qudt.org/vocab/unit/DeciM3",
                     "day": "https://qudt.org/vocab/unit/DAY",
                     "MPa": "http://qudt.org/vocab/unit/MegaPA",
                     "N/mm^2": "http://qudt.org/vocab/unit/N-PER-MilliM2",
                     "mm^2": "http://qudt.org/vocab/unit/MilliM2",
                     "g/mm^3": "http://qudt.org/vocab/unit/GM-PER-CentiM3",
                     "mm/mm": "http://qudt.org/vocab/unit/PER-MilliM",
                     "kg": "http://qudt.org/vocab/unit/KiloGM",
                     "MilliM": "http://qudt.org/vocab/unit/MilliM",
                     "GM": "http://qudt.org/vocab/unit/GM",
                     "KN": "https://qudt.org/vocab/unit/KN",
                     "KiloGM-PER-M3": "http://qudt.org/vocab/unit/KiloGM-PER-M3",
                     "KiloGM-PER-DeciM3": "http://qudt.org/vocab/unit/KiloGM-PER-DeciM3",
                     "DeciM3": "http://qudt.org/vocab/unit/DeciM3",
                     "DAY": "https://qudt.org/vocab/unit/DAY",
                     "MegaPA": "http://qudt.org/vocab/unit/MegaPA",
                     "N-PER-MilliM2": "http://qudt.org/vocab/unit/N-PER-MilliM2",
                     "MilliM2": "http://qudt.org/vocab/unit/MilliM2",
                     "GM-PER-CentiM3": "http://qudt.org/vocab/unit/GM-PER-CentiM3",
                     "PER-MilliM": "http://qudt.org/vocab/unit/PER-MilliM",
                     "KiloGM": "http://qudt.org/vocab/unit/KiloGM"}

    # Add more unit mappings as needed

    debug_counter = 0
    
    # Check if the unit string exists in the dictionary
    for key in input_metadata:
        if key.endswith("_Unit"):
            if input_metadata[key] in unit_mappings:
                input_metadata[key] = unit_mappings[input_metadata[key]]

            debug_counter += 1

    logger.debug("Replaced units:")
    logger.debug(str(debug_counter))

    output_metadata = input_metadata

    return output_metadata

def unit_conversion_json(input_metadata):
    """
    Script to convert the units extracted from the raw data into valid ontology
    instances.

    :param input_metadata: dictionary of metadata
    :return: dictionary of converted units
    """

    # Define the unit mappings as a dictionary
    unit_mappings = {"mm": "MilliM",
                     "g": "GM",
                     "kN": "KN",
                     "kg/m^3": "KiloGM-PER-M3",
                     "kg/dm^3": "KiloGM-PER-DeciM3",
                     "dm^3": "DeciM3",
                     "day": "DAY",
                     "MPa": "MegaPA",
                     "N/mm^2": "N-PER-MilliM2",
                     "mm^2": "MilliM2",
                     "g/mm^3": "GM-PER-CentiM3",
                     "mm/mm": "PER-MilliM",
                     "kg": "KiloGM"}

    # Add more unit mappings as needed

    debug_counter = 0
    
    # Check if the unit string exists in the dictionary
    for key in input_metadata:
        if key.endswith("Unit"):
            if input_metadata[key] in unit_mappings:
                input_metadata[key] = unit_mappings[input_metadata[key]]

            debug_counter += 1

    logger.debug("Replaced units:")
    logger.debug(str(debug_counter))

    output_metadata = input_metadata

    return output_metadata
