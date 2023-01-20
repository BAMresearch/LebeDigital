# Interbis

## Motivation

Interbis is an extension of the [Pybis](https://pypi.org/project/PyBIS/) library meant for improving data synchronisation with an openBIS datastore.

## Glossary

### Disclaimer:
The terms Sample and Objects are aliases for each other, you can use them synonymously within the pybis environment. We use Sample and Collection primarily in Interbis

- **Interbis** - A Subclass of the Openbis class in the pybis library. It is used for connecting to the openBIS Database and for all Database accesses.
- **Space** - The Highest "directory" abstraction in the Datastore. Contains all other openBIS "directories". One Space is defined per User or per Group.
- **Project** - Second highest "directory" abstraction in the Datastore. Contained in a Space, you can define infinitely many of them.
  - _Example_: UCT_Compression
- **Collection / Experiment** - Third highest "directory" abstraction in the Datastore. Contained in a Project, you can define infinitely many of them.
  - _Example_: Concrete_Compression_Testseries01
- **Sample / Object Type** - A blueprint defining the general structure of a specific experimental set-up (metadata, relations ... ), you need to define in order to upload Samples / Objects of real world experiments.
  - _Example_: EXPERIMENTAL_STEP_UCT
- **Sample / Object** - An abstract representation of your experiment with your defined metadata **also called Experimental Steps** in the openBIS Database. Contained in a Collection, you can define infinitely many of them.
  - _Example_: UCT_Compression_Experiment_17_10_2022
- **Dataset** - An openBIS abstraction containing your data files like Excel sheets or csvs with some metadata information.
  - _Example_: UCT_Compression_Experiment_17_10_2022_RAW_DATA
- **Property** - One piece of metadata you use do describe your Samples / Objects. 
  - _Example_: Name of your sample / object
- **Parents** of a Sample / Object - You can define another Sample / Object to be a parent of a different sample / object.
  - _Example_: A concrete mixture could be a **parent** of a concrete compression experiment
  - _Note_: Also defined are **Children** with the same relation in the other direction.

> You can look at the definition from pybis [here](https://pypi.org/project/PyBIS/)

## Structure

### Interbis

Is an extension of the existing Openbis class defined in pybis and as such extends all pybis functionality.
Interbis provides additional functionality for getting an overview of the datastore and find out what the structure of the Database looks like.
Interbis also makes it easier to search for existing samples / objects and retrieve their Identifiers from their metadata like Name.

## Example Workflow of Interbis

### 1. Connecting to the Database

```python
from lebedigital.openbis.interbis import Interbis

o = Interbis('https://yourdatabase.com') # Define an Interbis object with the url to your openBIS Database
o.connect_to_datastore() # The username will be read from your account (windows/linux) and you will have to input the password
```

### 2. Create a blueprint for the Samples / Objects - create a new **Sample / Object Type**

```python
# In order to define a new Sample / Object Type use the following method

o.create_sample_type(
    sample_code='EXPERIMENTAL_STEP_UCT', # Te identifier of the new sample / object type
    sample_prefix='UCL', # The prefix of the sample / object. Will appear before every sample / object code
    sample_properties={
        '$name': ['VARCHAR', 'Name', 'Name'] # Default system property
        'date': ['DATE', 'Date' 'Date of the experiment']
        'property1': ['VARCHAR', 'my_label1', 'my_description1'], # Define properties here with data in the order
        'property2': ['INTEGER', 'my_label2', 'my_description2'], # [Data Type, Label, Description]
        'property3': ['FLOAT', 'my_label2', 'my_description2']
    }
)
```
> You can find all the possible property types like 'VARCHAR' or 'FLOAT' [here](https://pypi.org/project/PyBIS/) under 'create property types'

> You can also define the Sample / Object Type in the Web-View of the Datastore. You need to be an admin to create new Sample / Object Types with both options.

### 3. Create and upload samples / objects from type "EXPERIMENTAL_STEP_UCT"

For each measured experiment a sample / object of the corresponding type (EXPERIMENTAL_STEP_UCT) has to be created.

#### a) Define the Sample

```python
concrete_experiment = o.new_sample(
    type = 'EXPERIMENTAL_STEP_UCT', # Has to be the same as a defined Sample Type. Here the one from chapter 2
    space = '<Your Username>', # your preassigned personal space
    project = 'Unilateral Conrete Compression',
    collection = '<Your Username>/Unilateral Conrete Compression/Concrete_Compression_Testseries01',
)

experiment_metadata = {'$name': 'Concrete_Experiment_1', 'date': '17:10:2022'}
concrete_experiment.set_props(experiment_metadata)
```

#### b) Upload the Sample to the Database
```python
sample.save()
```
#### Good to know:
* To upload a sample / object of a given type that type has to already exist in the datastore (See 2. Create a new Sample / Object Type).
* The name of the sample / object has to be the same as the `$name` parameter in the metadata dictionary in order to be able to find the sample / object efficiently again.
* You can look up which metadata the Experimental Step Type will accept by either getting the Excel import template running the Interbis `get_metadata_import_template()` method or get the properties of an Experimental Step Type by running the Interbis `get_sample_type_properties` method.

> Return of `o.get_metadata_import_template(EXPERIMENTAL_STEP_UCT)` as an Excel Sheet or Pandas Dataframe (o = connected Interbis instance)

|     | Param     | Label     | Description            | Value |
|-----|-----------|-----------|------------------------|-------|
| 0   | $name     | Name      | Name                   |       |
| 1   | date      | Date      | Date of the experiment |       |
| 2   | property1 | my_label1 | my_description1        |       |
| 3   | property2 | my_label2 | my_description2        |       |
| 4   | property3 | my_label3 | my_description3        |       |

### 4. Upload corresponding Datasets.

```python
concrete_dataset = o.new_dataset(
    type = 'RAW_DATA', # Or PROCESSED_DATA or some other Dataset Type. Check your own datastore with o.get_dataset_types() for available datasets
    collection = concrete_experiment.collection,
    sample = concrete_experiment.identifier,
    files = ['path_to_file1', 'path_to_file2'],
    props = {'$name': 'UCT_Compression_Experiment_17_10_2022_RAW_DATA',}
)
```

### If you have done everything correctly you should be able to see your Experimental Step of type 'EXPERIMENTAL_STEP_UCT' and your Datasets in the Web-View of the Database in the location you defined. 