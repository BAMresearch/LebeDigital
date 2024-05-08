## Here are the steps to set up the Flask project:

1. **Clone the Repository:** Clone the project repository from GitHub to your machine.
   ``` 
   git clone https://github.com/BAMresearch/LebeDigital.git
   ```

After cloning successfully, change to the target branch.
   ``` 
   git checkout workflowTest 
   ```

2. **Navigate to the Project Directory:** Go to the web folder inside the project directory.
    ```
    cd web
    ```

3. **Set Up a Virtual Environment:** It’s a good practice to create a virtual environment for your Python projects to isolate dependencies. You can use venv module which is included in standard Python library.
    ```
    python3 -m venv venv
    ```

4. **Activate the Virtual Environment:** Activate corresponding environment.

For Linux:
    ``` 
    source venv/bin/activate 
    ```
  
For Windows:
    ``` 
    .\\venv\\Scripts\\activate 
    ```

5. **Install Dependencies:** Install necessary dependencies using requirements.txt file.
    ``` 
    pip install -r requirements.txt
    ```

6. **Set Up Configuration:** Make copy of config.json.sample file, rename it config.json, and fill in necessary credentials and tokens.

Copy command:
   ``` 
   cp config.json.sample config.json 
   ```
   
   Open `config.json` and replace placeholders with actual credentials.

7. **Run Application:** Start Flask application which is defined in server.py file.
 
Run command:
     ``` 
     python server.py  
     ```

