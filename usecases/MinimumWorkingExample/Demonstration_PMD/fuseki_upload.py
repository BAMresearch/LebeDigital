import requests

def upload_ttl_file(fuseki_url, dataset_name, ttl_file_path):
    # Prepare the URL for uploading data to a specific dataset
    upload_url = f"{fuseki_url}/{dataset_name}/data"

    # Read the TTL file content
    with open(ttl_file_path, "rb") as file:
        ttl_data = file.read()

    # Define the headers for the request (content type should be "text/turtle")
    headers = {"Content-Type": "text/turtle"}

    # Send the POST request to upload the data
    response = requests.post(upload_url, data=ttl_data, headers=headers)

    if response.status_code == 200:
        print("TTL file upload successful.")
    else:
        print(f"Error occurred: {response.status_code} - {response.text}")

if __name__ == "__main__":
    # Replace the following with your Fuseki server's URL, dataset name, and TTL file path
    fuseki_url = "https://mechanics:7fewqg2nDWtZb5HC@fuseki.matolab.org"
    dataset_name = "testKG"
    ttl_file_path = "C:/develop/lebedigital-new/Lebedigital/usecases/MinimumWorkingExample/Mapping_Example/testMapped.ttl"

    upload_ttl_file(fuseki_url, dataset_name, ttl_file_path)
