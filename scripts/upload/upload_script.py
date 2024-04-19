import requests


def upload_binary_to_existing_docker(binary_data, config):
    """
    Uploads a ttl triple file in binary format to the dataset on the ontodocker,
    configured in the config file.

    :param binary_data: ttl triple file in binary format
    :param config: Config-File, containing configurations: DOCKER_TOKEN, ontodocker_url,
                    dataset_name and triplestore_server
    :return: if success 0, else 1
    """

    # Set Authorization Header with token
    headers = {"Authorization": f"Bearer {config['DOCKER_TOKEN']}", "Content-Type": "text/turtle"}

    # Prepare URL's for upload
    upload_url = f'{config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}/upload'
    print(f"Uploading data to dataset '{config['dataset_name']}' at '{upload_url}'.")

    # Upload the binary data to the dataset
    upload_response = requests.post(upload_url, headers=headers, data=binary_data)
    print(upload_response.content.decode())

    # Status based on response code
    if upload_response.status_code == 200:
        print("Data successfully uploaded.")
        return 0
    else:
        print("Failed to upload data.")
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
        results = response.json()
        return results
    else:
        print(f"Fehler bei der Anfrage: HTTP {response.status_code}, {response.text}")
        return f"Fehler bei der Anfrage: HTTP {response.status_code}"


def clear_dataset(config):
    """
    Clears the dataset selected in config

    :param config: Config-File, containing configurations: DOCKER_TOKEN, ontodocker_url,
                    dataset_name and triplestore_server
    :return: true
    """

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
    requests.post(f'{config["ontodocker_url"]}/api/{config["triplestore_server"]}/{config["dataset_name"]}/update',
                  headers=headers, params={"update": del_query}).content.decode()

    return
