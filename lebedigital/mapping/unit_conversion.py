# Script to convert the units extracted from the raw data into valid ontology 
# instances.


from loguru import logger


def unit_conversion(input_metadata):
    """

    Parameters:
    ----------
    input_metadata : dict
        Dictionary containing the extracted metadata.

    Returns:
    -------
    output_metadata : dict
        Dictionary with units replaced by link to an instance defined in the PMD core.
    """
    # Define the unit mappings as a dictionary
    unit_mappings = {
       'mm': 'http://qudt.org/vocab/unit/MilliM',
       'g': 'http://qudt.org/vocab/unit/GM',
       'kN': 'https://qudt.org/vocab/unit/KN',
       'kg/m³': 'http://qudt.org/vocab/unit/KiloGM-PER-M3',
       'kg/dm³': 'http://qudt.org/vocab/unit/KiloGM-PER-DeciM3',
       'dm³': 'http://qudt.org/vocab/unit/DeciM3'

        # Add more unit mappings as needed
    }


    debug_counter = 0
    # Check if the unit string exists in the dictionary
    for key in input_metadata:
        if key.endswith("_Unit"):
            if input_metadata[key] in unit_mappings:
                input_metadata[key] = unit_mappings[input_metadata[key]]  ### Replace the unit "mm" or something else to a proper ontology unit like "https://..."

            debug_counter += 1

    logger.debug("Replaced units:")
    logger.debug(str(debug_counter))

    output_metadata = input_metadata

    return output_metadata

