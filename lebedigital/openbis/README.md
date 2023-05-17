# Interbis

## Motivation

Interbis is an extension of the [Pybis](https://pypi.org/project/PyBIS/) library meant for improving data synchronisation with an openBIS datastore.

## Glossary

### Disclaimer:
The terms sample and objects are aliases for each other, you can use them synonymously within the pybis environment. We use Sample and Collection primarily in Interbis

- **Interbis** - A Subclass of the Openbis class in the pybis library. It is used for connecting to the openBIS database and for all database accesses.
- **Space** - The highest "directory" abstraction in the datastore. Contains all other openBIS "directories". One space is defined per user or per group.
- **Project** - Second highest "directory" abstraction in the datastore. Contained in a space, you can define infinitely many of them.
  - _Example_: UCT_Compression
- **Collection / Experiment** - Third highest "directory" abstraction in the datastore. Contained in a project, you can define infinitely many of them.
  - _Example_: UCT_Compression_Testseries01
- **Sample Type / Object Type** - A blueprint defining the general structure of a specific experimental set-up (metadata, relations ... ), you need to define in order to upload samples / objects of real world experiments.
  - _Example_: EXPERIMENTAL_STEP_UCT
- **Sample / Object** - An abstract representation of your experiment with your defined metadata **also called Experimental Step** in the openBIS database. Contained in a collection, you can define infinitely many of them.
  - _Example_: UCT_Compression_Experiment_17_10_2022
- **Dataset** - An openBIS abstraction containing your data files like Excel sheets, csvs or dats with some metadata information.
  - _Example_: UCT_Compression_Experiment_17_10_2022_RAW_DATA
- **Property** - One piece of metadata you use do describe your samples / objects. 
  - _Example_: NAME
- **Parents** of a sample / object - You can define another sample / object to be a parent of a different sample / object.
  - _Example_: A concrete mixture could be a **parent** of a concrete compression experiment
  - _Note_: Also defined are **Children** with the same relation in the other direction.

> You can look at the definition from pybis [here](https://pypi.org/project/PyBIS/)

### Indentification in openBIS

Each sample / object in openBIS will automatically ge
- a **Permid** - A unique number
- a **Code** - A unique name consists of the sample type's / object type's prefix plus a counting up number
- an **Identifier** - Path SPACE/PROJECT/CODE
- a **Path** - Directory path to current sample / object usually SPACE/PROJECT/COLLECTION/CODE.
- a **Type** - Name of the sample type / object type of current sample / object.
- a **Experiment** - Directory path to collection, if current sample / object is part of a collection / experiment.

The **permid** or the **code** can be used to clearly identify a specific sample / object and e.g. link it to others.

## Structure

### Interbis

Is an extension of the existing openbis class defined in pybis and as such extends all pybis functionality.
Interbis provides additional functionality for getting an overview of the datastore and find out what the structure of the database looks like.
Interbis also makes it easier to search for existing samples / objects and retrieve their identifiers from their metadata like NAME.

## Example Workflow of Interbis

### 1. Connecting to the database

```python
from lebedigital.openbis.interbis import Interbis

o = Interbis('https://yourdatabase.com') # Define an Interbis object with the url to your openBIS Database
o.connect_to_datastore(username='name', password=None) # If username is not given, the username will be read from your account (windows/linux). If the password is not given, there will be a prompt.
```

### 2. Create a blueprint for the samples / objects - create a new **sample type / object type**

```python
# In order to define a new sample / object type use the following method (only possible as admin)

o.create_sample_type(
    sample_code='EXPERIMENTAL_STEP_UCT', # The code identifier of the new sample / object type
    sample_prefix='UCT', # The prefix of the sample / object. Will appear before every sample / object code
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

> You can also define the sample / object type in the Web-View of the datastore. 

> You need to be an admin to create new sample / object types.

### 3. Create and upload samples / objects from type "EXPERIMENTAL_STEP_UCT"

For each measured experiment a sample / object of the corresponding type (EXPERIMENTAL_STEP_UCT) has to be created.

#### a) Define the sample

```python
concrete_experiment = o.new_sample(
    type = 'EXPERIMENTAL_STEP_UCT', # Has to be the same as a defined sample type. Here the one from chapter 2
    space = '<Your Username>', # your preassigned personal space
    project = 'UCT_Compression', 
    collection = '<Your Username>/UCT_Compression/UCT_Compression_Testseries01',
)

experiment_metadata = {'$name': 'Concrete_Experiment_1', 'date': '17:10:2022'}
concrete_experiment.set_props(experiment_metadata)
```

#### b) Upload the Sample to the Database
```python
sample.save()
```
#### Good to know:
* To upload a sample / object of a given type that type has to already exist in the datastore (See 2. Create a new sample type / object type).
* The project and collection the sample should be uploaded to, must exist before uploading. Creating new projects/ collection can be done via the Web-View (ELN) or by pybis:
```python
project_obj = o.new_project(
                            space='<Your Username>', 
                            code='UCT_Compression',
                            description="...")
project_obj.save()
collection_obj = o.new_collection(
                            project='UCT_Compression', 
                            code='UCT_Compression_Testseries01', 
                            type="COLLECTION")
collection_obj.save()
```
* You can look up which metadata the Experimental Step Type will accept by either getting the Excel import template running the Interbis `get_metadata_import_template()` method or get the properties of an Experimental Step Type by running the Interbis `get_sample_type_properties` method. (Or via the admin Web-View.)

> Return of `o.get_metadata_import_template(EXPERIMENTAL_STEP_UCT)` as an Excel Sheet or Pandas Dataframe (o = connected Interbis instance)

|     | Param     | Label     | Description            | Value |
|-----|-----------|-----------|------------------------|-------|
| 0   | $name     | Name      | Name                   |       |
| 1   | date      | Date      | Date of the experiment |       |
| 2   | property1 | my_label1 | my_description1        |       |
| 3   | property2 | my_label2 | my_description2        |       |
| 4   | property3 | my_label3 | my_description3        |       |

> If you use the Excel Sheet above for filling in metadata, then the upload to the sample can be done via `o.import_props_from_template(<path to filled out Excell sheet>, sample)`.

### 4. Upload corresponding datasets.

```python
concrete_dataset = o.new_dataset(
    type = 'RAW_DATA', # Or PROCESSED_DATA or some other Dataset Type. Check your own datastore with o.get_dataset_types() for available datasets
    collection = concrete_experiment.collection,
    sample = concrete_experiment.identifier,
    files = ['path_to_file1', 'path_to_file2'],
    props = {'$name': 'UCT_Compression_Experiment_17_10_2022_RAW_DATA',}
)
```

### If you have done everything correctly you should be able to see your Experimental Step of type 'EXPERIMENTAL_STEP_UCT' and your datasets in the Web-View of the database in the location you defined. 

### 5. Searching for specific samples / objects in database.

```python
# a Pandas DataFrame of all suitable samples / objects
samples = o.get_samples(
    space="<Your Username>", # search in this space
    project="TEST_AMBETON", # search in this project 
    props=['$NAME', 'DATE'], # show these properties in the results using props="*" retrieves all properties
    where={
        "DATE": "2022-05-10", # query
    },
).df
```
> For detailed information see [here](https://pypi.org/project/PyBIS/) under 'search for samples / objects'.
> Searching is also possible via Web-View.

## Generating and working with type checkers.

> You can generate type-checkers based on pydantic models to check your samples for formatting errors before uploading them to openbis

```python
sample_properties = {"$name": "sample_name"}

SampleModel = o.generate_typechecker("SAMPLE_TYPE")
sample_model_return = SampleModel(**sample_properties)
sample_properties = sample_model_return.dict(exclude_unset=True) # exclude_unset is for not including other property types that were not set in sample_properties
```

> The sample_properties have been checked for their types and casted/formatted if possible. The functionality includes:

- Casting strings like "2.13" or "7" to float or integer respectively
- Checking if CONTROLLED_VOCABULARY Properties are withing their defined vocabularies
- Casting date-strings into the required format. Ex "17.10.2023 10:45" will be formatted to "2023-10-17"
