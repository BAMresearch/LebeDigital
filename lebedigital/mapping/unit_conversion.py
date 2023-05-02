# Script to convert the units extracted from the raw data into valid ontology 
# instances.

from loguru import logger


def unit_conversion(input_metadata):
    """
    Converts the units given in a metadata-file into ontology instances for
    further mapping process.

    Parameters:
    ----------
    input_metadata : dict
        Dictionary containing the extracted metadata.

    Returns:
    -------
    output_metadata : dict
        Dictionary with units replaced by link to an instance defined in in the PMD core.
    """

    keys = list(input_metadata.keys()) 
    debug_counter = 0

    for key in keys:
        if key[-4:] == "Unit":
            input_metadata[key] = "Successfully replaced this unit!" #placeholder for testing the conversion until I get a table with real data
            debug_counter += 1

    logger.debug("Replaced units:")
    logger.debug(str(debug_counter))

    output_metadata = input_metadata

    return output_metadata