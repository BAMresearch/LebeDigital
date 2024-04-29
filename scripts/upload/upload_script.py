import requests
from loguru import logger


def upload_binary_to_existing_docker(binary_data, config):
    """
    Uploads a ttl triple file in binary format to the dataset on the ontodocker,
    configured in the config file.

    :param binary_data: ttl triple file in binary format
    :param config: Config-File, containing configurations: DOCKER_TOKEN, ontodocker_url,
                    dataset_name and triplestore_server
    :return: if success 0, else 1
    """

    logger.debug("-" * 20)
    logger.debug(f"Startet uploading to: {config['ontodocker_url']} for file with length: {len(binary_data)}")
    logger.debug(
        f'Configuration: {config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}')

    # Set Authorization Header with token
    headers = {"Authorization": f"Bearer {config['DOCKER_TOKEN']}", "Content-Type": "text/turtle"}

    # Prepare URL's for upload
    upload_url = f'{config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}/upload'

    # Upload the binary data to the dataset
    upload_response = requests.post(upload_url, headers=headers, data=binary_data)
    logger.debug(upload_response.content.decode())

    # Status based on response code
    if upload_response.status_code == 200:
        logger.debug("Data successfully uploaded.")
        logger.debug("-" * 20)
        return 0
    else:
        logger.error("Failed to upload data to Ontodocker.")
        logger.debug("-" * 20)
        return 1


def send_sparql_query(query, config):
    """
    Sends a SPARQL query to a Jena Fuseki server (Ontodocker at BAM) and checks
    whether the search was successful. If so, the search result is transferred.

    :param query: The SPARQL-Query as String
    :param config: Config-File, containing configurations: DOCKER_TOKEN, ontodocker_url,
                    dataset_name and triplestore_server
    :return: result of Query if successful, else error
    """

    logger.debug("-" * 20)
    logger.debug(f'Sending Sparql-Query to Ontodocker: {query}')
    logger.debug(f'Configuration: {config["ontodocker_url"]}/api/{config["triplestore_server"]}/'
                 f'{config["dataset_name"]}'
                 f'/sparql')

    # Set Authorization Header with token
    headers = {
        "Authorization": f"Bearer {config['DOCKER_TOKEN']}"
    }

    # Add the Query to data
    data = {
        'query': query
    }

    # Send request to the server
    response = requests.get(f'{config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}'
                            f'/sparql',
                            headers=headers, params=data)

    # Return based on response code
    if response.status_code == 200:
        logger.debug("Query successful.")
        logger.debug("-" * 20)
        results = response.json()
        return results
    else:
        logger.error(f"Error in the query request for the Ontodocker: HTTP {response.status_code}, {response.text}")
        return f"Fehler bei der Anfrage: HTTP {response.status_code}"


def clear_dataset(config):
    """
    Clears the dataset selected in config

    :param config: Config-File, containing configurations: DOCKER_TOKEN, ontodocker_url,
                    dataset_name and triplestore_server
    :return: true
    """

    logger.debug("-" * 20)
    logger.debug("Clearing Dataset and Ontodocker.")
    logger.debug(f'Configuration: {config["ontodocker_url"]}/api/{config["triplestore_server"]}/'
                 f'{config["dataset_name"]}/update')

    # Set Authorization Header with token
    headers = {
        "Authorization": f"Bearer {config['DOCKER_TOKEN']}"
    }

    # Set the deletion Query
    del_query = """
        DELETE WHERE
        {
          ?s ?p ?o .
        }
        """

    # Post the deletion request to the server
    requests.post(f'{config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}'
                  f'/update',
                  headers=headers, params={"update": del_query}).content.decode()

    logger.debug("Deletion request successfully sent.")
    logger.debug("-" * 20)

    return
