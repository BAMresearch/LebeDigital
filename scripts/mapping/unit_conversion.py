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
                     "GPa": "https://qudt.org/vocab/unit/GigaPA",
                     "cylindrical": "https://w3id.org/cpto/Cylinder",
                     "cubical":  "https://w3id.org/cpto/Cube",
                     "N/mm^2": "http://qudt.org/vocab/unit/N-PER-MilliM2",
                     "mm^2": "http://qudt.org/vocab/unit/MilliM2",
                     "g/mm^3": "http://qudt.org/vocab/unit/GM-PER-CentiM3",
                     "mm/mm": "http://qudt.org/vocab/unit/PER-MilliM",
                     "kg": "http://qudt.org/vocab/unit/KiloGM"}

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
