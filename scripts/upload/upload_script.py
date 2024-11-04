from loguru import logger
from rdflib import Graph
import requests
from requests.auth import HTTPBasicAuth
import os
from urllib.parse import urljoin

# Scripts for handling the ttl dataset

def check_if_file_exists(config):

    try:
        with open(config["dataset_name"], "r") as file:
            return True
    except FileNotFoundError:
        # create empty file
        with open(config["dataset_name"], "w") as file:
            file.write("")


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
        # Construct the correct upload URL for Fuseki
        base_url = config["fuseki"]["endpoint_url"].rstrip('/')  # Remove trailing slash if present
        dataset = config["fuseki"]["dataset_name"]
        
        # Use the correct endpoint structure as observed in network request
        upload_url = f"{base_url}/{dataset}/data"
        
        logger.debug(f"Attempting to upload to URL: {upload_url}")
        
        # Set up authentication
        auth = HTTPBasicAuth(
            config["fuseki"]["username"],
            config["fuseki"]["password"]
        )

        # Set up headers
        headers = {
            'Content-Type': 'text/turtle'
        }

        # Make the POST request as observed in network request
        response = requests.post(
            upload_url,
            data=binary_data,
            headers=headers,
            auth=auth,
            verify=True
        )

        # Check if the upload was successful
        if response.status_code == 200:
            logger.info(f"Successfully uploaded TTL data to Fuseki server")
            return True
        else:
            logger.error(f"Failed to upload TTL data. Status code: {response.status_code}")
            logger.error(f"Response: {response.text}")
            return False

    except Exception as e:
        logger.error(f"Error uploading TTL data to Fuseki: {str(e)}")
        return False
    

def send_sparql_query(query, config):
    """
    Query the Fuseki dataset with a SPARQL-Query

    :param query: The SPARQL-Query as String
    :param config: Config-File, containing Fuseki settings
    :return: result of Query if successful, else error
    """
    logger.debug("-" * 20)
    logger.debug(f'Sending Sparql-Query: {query}')

    try:
        # Construct the SPARQL endpoint URL
        base_url = config["fuseki"]["endpoint_url"].rstrip('/')
        dataset = config["fuseki"]["dataset_name"]
        query_url = f"{base_url}/{dataset}/sparql"

        # Set up authentication
        auth = HTTPBasicAuth(
            config["fuseki"]["username"],
            config["fuseki"]["password"]
        )

        # Set up headers
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
            'Accept': 'application/json'
        }

        # Set up parameters
        params = {
            'query': query
        }

        # Make the POST request
        response = requests.post(
            query_url,
            headers=headers,
            auth=auth,
            data=params,
            verify=True
        )

        if response.status_code == 200:
            # Parse JSON response
            json_response = response.json()
            
            # Extract variables (column names) from the response
            if 'head' in json_response and 'vars' in json_response['head']:
                column_order = json_response['head']['vars']
            else:
                column_order = []

            # Create results dictionary with the same structure as before
            results_dict = {
                'columns': column_order,
                'data': []
            }

            # Process bindings (results)
            if 'results' in json_response and 'bindings' in json_response['results']:
                for binding in json_response['results']['bindings']:
                    row_dict = {}
                    for col in column_order:
                        # Check if the variable exists in the binding
                        if col in binding:
                            row_dict[col] = binding[col]['value']
                        else:
                            row_dict[col] = ''
                    results_dict['data'].append(row_dict)

            logger.debug("Query successful.")
            logger.debug(results_dict)
            return results_dict
        else:
            logger.error(f"Query failed with status code: {response.status_code}")
            logger.error(f"Response: {response.text}")
            return None

    except Exception as e:
        logger.error(f"Error in the query request: {str(e)}")
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
