from loguru import logger
from rdflib import Graph

# Scripts for handling the ttl dataset

def check_if_file_exists(config):

    try:
        with open(config["dataset_name"], "r") as file:
            return True
    except FileNotFoundError:
        # create empty file
        with open(config["dataset_name"], "w") as file:
            file.write("")


def upload_binary_to_existing_docker(binary_data, config):
    """
    Adds a ttl triple file in binary format to the dataset "datastore.ttl"

    :param binary_data: ttl triple file in binary format
    :return: if success 0, else 1
    """

    check_if_file_exists(config)

    logger.debug("-" * 20)

    try:
        # Create a Fraph for the current dataset
        dataset_graph = Graph()
        
        # Load the current dataset from the file
        dataset_graph.parse(config["dataset_name"], format="turtle")
        
        # Create a new Graph for the new data
        new_graph = Graph()
        
        # Load the new data from the binary file
        new_graph.parse(data=binary_data, format="turtle")
        
        # Add the new data to the current dataset
        dataset_graph += new_graph
        
        # Save the updated dataset to the file
        dataset_graph.serialize(destination=config["dataset_name"], format="turtle")
        
        logger.debug("Data successfully added to the dataset.")

        return 0
    
    except Exception as e:
        logger.error(f"Failed to add data to the dataset: {e}")

        return 1


def send_sparql_query(query, config):
    """
    Query the dataset with a SPARQL-Query

    :param query: The SPARQL-Query as String
    :param config: Config-File, containing the dataset name
    :return: result of Query if successful, else error
    """

    check_if_file_exists(config)

    logger.debug("-" * 20)
    logger.debug(f'Sending Sparql-Query: {query}')

    try:
        # Erstellen eines Graphen für das bestehende Dataset
        dataset_graph = Graph()
        
        # Laden des bestehenden Datasets
        dataset_graph.parse(config["dataset_name"], format="turtle")
        
        # Ausführen der SPARQL-Abfrage
        results = dataset_graph.query(query)
        
        # Ergebnisse in einer Liste speichern
        # Konvertieren Sie die Ergebnisse in ein Dictionary
        results_dict = []
        for row in results:
            row_dict = {str(var): str(row[var]) for var in results.vars}
            results_dict.append(row_dict)
        
        logger.debug("Query successful.")
        logger.debug(results_dict)
        # Rückgabe der Ergebnisse
        return results_dict
    
    except Exception as e:
        logger.error(f"Error in the query request: {e}")

        return None


def clear_dataset(config):
    """
    Clears the dataset selected in config

    :param config: Config-File, containing configurations: DOCKER_TOKEN, ontodocker_url,
                    dataset_name and triplestore_server
    :return: true
    """

    check_if_file_exists(config)

    logger.debug("-" * 20)
    logger.debug("Clearing Dataset and Ontodocker.")

    try:
        with open (config["dataset_name"], "w") as file:
            file.write("")
        
        logger.debug("Deletion successfull.")
        logger.debug("-" * 20)
    except Exception as e:
        logger.error(f"Error in clearing the dataset: {e}")

    return
