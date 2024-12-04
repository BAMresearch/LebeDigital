from loguru import logger
import requests
from requests.auth import HTTPBasicAuth
from io import BytesIO

# Constants for HTTP status codes and content types
HTTP_SUCCESS_CODES = [200, 201, 204]
CONTENT_TYPE_TURTLE = 'text/turtle'
CONTENT_TYPE_FORM = 'application/x-www-form-urlencoded'
CONTENT_TYPE_JSON = 'application/json'

def create_auth(config):
    """
    Create HTTP Basic Authentication object from config
    
    Args:
        config: Configuration dictionary containing Fuseki credentials
    Returns:
        HTTPBasicAuth object
    """
    return HTTPBasicAuth(
        config["fuseki"]["username"],
        config["fuseki"]["password"]
    )


def construct_fuseki_url(config, endpoint_type):
    """
    Construct Fuseki endpoint URL based on type
    
    Args:
        config: Configuration dictionary containing Fuseki settings
        endpoint_type: Type of endpoint ('data', 'sparql', 'update')
    Returns:
        str: Constructed URL
    """
    base_url = config["fuseki"]["endpoint_url"].rstrip('/')
    dataset = config["fuseki"]["dataset_name"]
    return f"{base_url}/{dataset}/{endpoint_type}"


def upload_ttl_to_fuseki(binary_data, config):
    """
    Upload TTL data to a Fuseki server with authentication
    
    Args:
        binary_data: The TTL file content
        config: Configuration dictionary containing Fuseki settings
    Returns:
        bool: True if upload was successful, False otherwise
    """
    try:
        upload_url = construct_fuseki_url(config, 'data')
        logger.debug(f"Attempting to upload to URL: {upload_url}")
        
        headers = {'Content-Type': CONTENT_TYPE_TURTLE}
        
        response = requests.post(
            upload_url,
            data=binary_data,
            headers=headers,
            auth=create_auth(config),
            verify=True
        )

        if response.status_code in HTTP_SUCCESS_CODES:
            logger.info("Successfully uploaded TTL data to Fuseki server")
            return True
        
        logger.error(f"Failed to upload TTL data. Status code: {response.status_code}")
        logger.error(f"Response: {response.text}")
        return False

    except Exception as e:
        logger.error(f"Error uploading TTL data to Fuseki: {str(e)}")
        return False


def send_sparql_query(query, config):
    """
    Query the Fuseki dataset with a SPARQL-Query

    Args:
        query: The SPARQL-Query as String
        config: Config-File containing Fuseki settings
    Returns:
        dict: Query results if successful, None if error
    """
    logger.debug("-" * 20)
    logger.debug(f'Sending Sparql-Query: {query}')

    try:
        query_url = construct_fuseki_url(config, 'sparql')
        
        headers = {
            'Content-Type': CONTENT_TYPE_FORM,
            'Accept': CONTENT_TYPE_JSON
        }

        response = requests.post(
            query_url,
            headers=headers,
            auth=create_auth(config),
            data={'query': query},
            verify=True
        )

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
        json_response: JSON response from Fuseki server
    Returns:
        dict: Processed results with columns and data
    """
    # Extract variables (column names) from the response
    column_order = json_response.get('head', {}).get('vars', [])
    
    results_dict = {
        'columns': column_order,
        'data': []
    }

    # Process bindings (results)
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
    Clears all triples from the Fuseki dataset
    
    Args:
        config: Config-File containing Fuseki settings
    Returns:
        bool: True if successful, False otherwise
    """
    logger.debug("-" * 20)
    logger.debug("Clearing Fuseki Dataset")

    try:
        update_url = construct_fuseki_url(config, 'update')
        logger.debug(f"Update endpoint: {update_url}")

        headers = {'Content-Type': CONTENT_TYPE_FORM}

        clear_query = """
        DELETE {
            ?s ?p ?o
        }
        WHERE {
            ?s ?p ?o
        }
        """

        response = requests.post(
            update_url,
            headers=headers,
            auth=create_auth(config),
            data={'update': clear_query},
            verify=True
        )

        if response.status_code in HTTP_SUCCESS_CODES:
            logger.debug("Dataset cleared successfully")
            return True
        
        logger.error(f"Failed to clear dataset. Status code: {response.status_code}")
        logger.error(f"Response: {response.text}")
        return False

    except Exception as e:
        logger.error(f"Error clearing dataset: {str(e)}")
        return False

def get_backup(config):
    """
    Get a Backup of the database

    Returns:
        ttl of database
    """
    logger.debug("-" * 20)
    logger.debug(f'Getting backup:')

    try:
        query_url = construct_fuseki_url(config, 'get')
        
        headers = {
            'Content-Type': CONTENT_TYPE_FORM,
            'Accept': CONTENT_TYPE_JSON
        }

        response = requests.get(
            query_url,
            headers=headers,
            auth=create_auth(config),
            verify=True
        )

        if response.status_code == 200:
            backup_stream = BytesIO(response.content)
            logger.debug("Backup erfolgreich abgerufen.")
            return backup_stream
        
        logger.error(f"Query failed with status code: {response.status_code}")
        logger.error(f"Response: {response.text}")
        return None

    except Exception as e:
        logger.error(f"Error in the backup request: {str(e)}")
        return None
