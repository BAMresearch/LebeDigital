import json
import requests
import loguru


def upload_binary_to_existing_docker(token, dataset_name, binary_data, triplestore="jena"):
    ontodocker_url = 'https://ontodocker-pmd.bam.de'  # Setze hier deine spezifische OntoDocker URL
    ontodocker_jwt = token

    # Setze den Authorization Header mit dem JWT-Token
    headers = {"Authorization": f"Bearer {ontodocker_jwt}", "Content-Type": "text/turtle"}

    # URL vorbereiten für den Upload
    upload_url = f'{ontodocker_url}/api/{triplestore}/{dataset_name}/upload'
    print(f"Uploading data to dataset '{dataset_name}' at '{upload_url}'.")

    # Lade die binäre Datei in das Dataset hoch
    upload_response = requests.post(upload_url, headers=headers, data=binary_data)
    print(upload_response.content.decode())

    # Statusmeldung basierend auf dem Antwortstatus
    if upload_response.status_code == 200:
        print("Data successfully uploaded.")
        return 1
    else:
        print("Failed to upload data.")
        return 0


def upload_to_docker(token, datasets):

    ontodocker_url = 'https://ontodocker-pmd.bam.de'  # replace with url to own ontodocker url.
    ontodocker_jwt = token

    # set ontodocker_jwt for further communication
    headers = {"Authorization": f"Bearer {ontodocker_jwt}"}

    # get all dataset ids
    print(requests.get(f'{ontodocker_url}/api/ds/all', headers=headers).content.decode())
    #exit()

    # create dataset ("jena" can be replaced with "blazegraph")
    dataset_name = "Test" #rename as desired
    triplestore = "jena"  #either "jena" or "blazegraph"
    print(requests.put(f'{ontodocker_url}/api/{triplestore}/{dataset_name}', headers=headers).content.decode())


    # destroy dataset
    #print(requests.delete(f'{ontodocker_url}/api/{triplestore}/{dataset_name}', headers=headers).content.decode())

    for entry in datasets:
        # upload file to dataset
        ontology_path = entry
        with open(ontology_path, 'rb') as data:
            headers.update({"Content-Type": "text/turtle"})
            print("url: ", f'{ontodocker_url}/api/{triplestore}/{dataset_name}/upload')
            print("data: ", data)
            print(requests.post(f'{ontodocker_url}/api/{triplestore}/{dataset_name}/upload', headers=headers, data=data).content.decode())

    # select
    #select_query = "SELECT ?s ?p ?o WHERE {?s ?p ?o}"
    #print(requests.get(f'{ontodocker_url}/api/{triplestore}/{dataset_name}/query', headers=headers, params={"query": select_query}).content.decode())

    # insert
    #insert_query = """
    #    PREFIX ns3: <http://qudt.org/schema/qudt/>
    #    INSERT DATA
    #    {
    #      <http://example/book1> dc:title "A new book" ;
    #                             dc:creator "A.N.Other" .
    #    }
    #    """
    #print(requests.post(f'{ontodocker_url}/api/{triplestore}/{dataset_name}/update', headers=headers, params={"update": insert_query}).content.decode())

    # delete
    #delete_query = """
    #    PREFIX dc: <http://purl.org/dc/elements/1.1/>
    #    DELETE DATA
    #    {
    #      <http://example/book1> dc:title "A new book" ;
    #                             dc:creator "A.N.Other" .
    #    }
    #    """
    #print(requests.post(f'{ontodocker_url}/api/{triplestore}/{dataset_name}/update', headers=headers, params={"update": delete_query}).content.decode())


    ## via SPARQLWrapper
    #from SPARQLWrapper import SPARQLWrapper, JSON
    #sparql = SPARQLWrapper(f"{ontodocker_url}/api/{triplestore}/{dataset_name}/sparql")

    # POST
    #insert_query = """
    #        PREFIX co: <https://w3id.org/pmd/co/>>
    #        INSERT DATA
    #        {
    #          <http://example/book1> dc:title "A new book" ;
    #                                 dc:creator "A.N.Other" .
    #        }
    #        """
    #sparql.setMethod("POST")
    #sparql.setQuery(insert_query)
    #sparql.addCustomHttpHeader("Authorization", f"Bearer {ontodocker_jwt}")
    #results = sparql.query().convert()
    #print(results)

    # GET
    #select_query = "SELECT ?s ?p ?o WHERE {?s ?p ?o}"
    #sparql.setMethod("GET")
    #sparql.setQuery(select_query)
    #sparql.addCustomHttpHeader("Authorization", f"Bearer {ontodocker_jwt}")
    #sparql.setReturnFormat(JSON)
    #results = sparql.query().convert()
    #for result in results["results"]["bindings"]:
    #   print(result["s"], result["p"], result["o"])


def send_sparql_query(query, token):
    """
    Sendet eine SPARQL-Query an einen Jena Fuseki Server und überprüft,
    ob die Suche erfolgreich war.

    :param endpoint_url: URL des SPARQL-Endpoints
    :param query: Die SPARQL-Query als String
    :return: True, wenn die Suche erfolgreich war, sonst False
    """

    ontodocker_url = 'https://ontodocker-pmd.bam.de'  # replace with url to own ontodocker url.
    ontodocker_jwt = token

    headers = {
        "Authorization": f"Bearer {ontodocker_jwt}"
    }

    data = {
        'query': query
    }

    dataset_name = "Test"
    triplestore = "jena"
    response = requests.get(f'{ontodocker_url}/api/{triplestore}/{dataset_name}/sparql', headers=headers, params=data)

    if response.status_code == 200:
        results = response.json()
        return results
    else:
        print(f"Fehler bei der Anfrage: HTTP {response.status_code}, {response.text}")
        return f"Fehler bei der Anfrage: HTTP {response.status_code}"
