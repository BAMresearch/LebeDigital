import requests
from loguru import logger
from SPARQLWrapper import SPARQLWrapper, POST

from rdflib import Graph, URIRef

def construct_ontodocker_url(config, endpoint_type):
    """
    Construct Ontodocker endpoint URL based on type
    
    Args:
        config: Configuration dictionary containing Ontodocker settings
        endpoint_type: Type of endpoint ('data', 'sparql', 'update')
    Returns:
        str: Constructed URL
    """
    base_url = config["ontodocker_url"].rstrip('/')
    return f"{base_url}/api/{config['triplestore_server']}/{config['dataset_name']}/{endpoint_type}"

def create_headers(config, content_type=None):
    """
    Create headers with authorization token
    
    Args:
        config: Configuration dictionary containing API key
        content_type: Optional content type to include in headers
    Returns:
        dict: Headers dictionary
    """
    headers = {"Authorization": f"Bearer {config['DOCKER_TOKEN']}"}
    if content_type:
        headers['Content-Type'] = content_type
    return headers


def extract_specific_triples_from_ttl(binary_data):
    """
    Extracts specific triples from a ttl triple file in binary format.

    :param binary_data: ttl triple file in binary format
    :param triples_to_extract: List of triples to extract
    :return: List of extracted triples
    """

    # Parse the ttl file
    g = Graph()
    g.parse(data=binary_data, format="turtle")

    # extract all triples from the ttl file
    triples_to_extract = [(str(s), str(p), str(o)) for s, p, o in g]

    return triples_to_extract

def delete_specific_triples_from_endpoint(ttl_binary, config):
    """
    Deletes specific triples from a SPARQL endpoint.

    :param ttl_binary: ttl triple file in binary format
    :param config: Config-File, containing configurations: SPARQL_ENDPOINT, DOCKER_TOKEN
    :return: if success 0, else 1
    """
    # Set Authorization Header with token
    headers = {
        "Authorization": f"Bearer {config['DOCKER_TOKEN']}"
    }

    triples_to_delete = extract_specific_triples_from_ttl(ttl_binary)

    # Prepare the SPARQL DELETE DATA query
    triples_str = "\n".join(f"<{s}> <{p}> <{o}> ." for s, p, o in triples_to_delete)
    query = f"""
        DELETE DATA {{
            GRAPH <{config['dataset_name']}> {{
                {triples_str}
            }}
        }}
    """

    logger.info(query)
    # Post the deletion request to the server
    response = requests.post(f'{config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}'
                  f'/update',
                  headers=headers, params={"update": query}).content.decode()

    # Send the SPARQL DELETE DATA query to the endpoint
    #sparql = SPARQLWrapper(config['SPARQL_ENDPOINT'])
    #sparql.setMethod(POST)
    #sparql.setQuery(query)
    #sparql.addParameter("Authorization", f"Bearer {config['DOCKER_TOKEN']}")
    logger.debug(response)

    return 0

def upload_ttl_to_ontodocker(binary_data, config):
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
    Query the Ontodocker dataset with a SPARQL-Query

    Args:
        query: The SPARQL-Query as String
        config: Config-File containing Ontodocker settings
    Returns:
        dict: Query results if successful, None if error
    """
    logger.debug("-" * 20)
    logger.debug(f'Sending Sparql-Query: {query}')

    try:
        query_url = construct_ontodocker_url(config, 'sparql')
        headers = create_headers(config)
        
        response = requests.get(
            query_url,
            headers=headers,
            params={'query': query},
            verify=True
        )

        print(response.json)
        if response.status_code == 200:
            return process_query_results(response.json())
        
        logger.error(f"Query failed with status code: {response.status_code}")
        logger.error(f"Response: {response.text}")
        return None

    except Exception as e:
        logger.error(f"Error in the query request: {str(e)}")
        return None

def process_query_results(json_response):
    """
    Process SPARQL query JSON response into structured format
    
    Args:
        json_response: JSON response from Ontodocker server
    Returns:
        dict: Processed results with columns and data
    """
    column_order = json_response.get('head', {}).get('vars', [])
    
    results_dict = {
        'columns': column_order,
        'data': []
    }

    bindings = json_response.get('results', {}).get('bindings', [])
    for binding in bindings:
        row_dict = {}
        for col in column_order:
            row_dict[col] = binding.get(col, {}).get('value', '')
        results_dict['data'].append(row_dict)

    logger.debug("Query successful.")
    logger.debug(results_dict)
    return results_dict


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

    return True
