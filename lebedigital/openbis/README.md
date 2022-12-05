# Interbis/ExpStep

## Motivation

### Interbis/Expstep is an extension of the Pybis library meant for improving data synchronisation with an openBIS datastore.

## Glossary

- **Interbis** - A Subclass of the Openbis class in the pybis library. It is used for connecting to the openBIS Database and for all Database queries like getting an overview of the Database
- **ExpStep** - A Class for working with openBIS Experimental Steps. Used for defining Objects which describe real Experiments that you want to describe with metadata and upload to the openBIS Database
- **Space** - The Highest "directory" abstraction in the Datastore. Contains all other openBIS "directories". One Space is defined per User or per Group.
  - _Example_: Some Paper about compressing concrete
- **Project** - Second highest "directory" abstraction in the Datastore. Contained in a Space, you can define infinitely many of them.
  - _Example_: Concrete_Compression
- **Collection / Experiment** - Third highest "directory" abstraction in the Datastore. Contained in a Project, you can define infinitely many of them.
  - _Example_: Concrete_Compression_Tests
- **Sample / Object** - An abstract representation of your experiment with your defined metadata **also called Experimental Steps** in the openBIS Database. Contained in a Collection, you can define infinitely many of them.
  - _Example_: Concrete_Compression_Experiment_17_10_2022
- **Dataset** - An openBIS abstraction containing your data files like Excel sheets or csvs.
  - _Example_: Concrete_Compression_Experiment_17_10_2022_RAW_DATA
- **Sample / Object Type** - A blueprint you need to define in order to upload Samples/Objects of real world experiments.
  - _Example_: EXPERIMENTAL_STEP_CONCRETE
- **Property** - One piece of metadata you use do describe your Samples/Objects. 
  - _Example_: Name of your sample
- **Parents** of a Sample - You can define another Sample to be a parent of a different sample.
  - _Example_: A concrete mixture could be a **parent** of a concrete compression experiment

> You can look at the definition from pybis [here](https://pypi.org/project/PyBIS/)

## Structure

### Interbis

Is an extension of the existing Openbis class defined in pybis and as such extends all pybis functionality.
Interbis provides additional functionality for getting an overview of the datastore and find out what the structure of the Database looks like.
Interbis also makes it easier to search for existing samples and retrieve their Identifiers from their metadata like Name.

### ExpStep

Is a Class which provides functionality to define Samples as they are accepted by pybis and in turn the openBIS datastore.
Example uses of the Class include uploading and retrieving samples and datasets from the Database.
I recommend deriving the class in your personal workflow and tweaking it to make it easier to write data processing and uploading scripts. 


## Example Workflow of Interbis/ExpStep - Concrete Measuring Project

### Start with connecting yourself to the Database

```python
from lebedigital.openbis.interbis import Interbis

o = Interbis('https://yourdatabase.com')
o.connect_to_datastore()
```

### Then you can define a blueprint from which you will create the Experimental Steps you have done in your research - create a new **Sample Type**

```python
# In order to define a new Sample Type use the following method

o.create_sample_type(
    sample_code='EXPERIMENTAL_STEP_CONCRETE_TEST',
    sample_prefix='CONCRETE_TEST',
    sample_properties={
        'property1': ['VARCHAR', 'my_label1', 'my_description1'],
        'property2': ['INTEGER', 'my_label2', 'my_description2'],
        'property3': ['FLOAT', 'my_label2', 'my_description2']
    }
)
```
> You can find all the possible property types like 'VARCHAR' or 'FLOAT' [here](https://pypi.org/project/PyBIS/) under 'create property types'

> You can also define the Sample Type in the Web-View of the Datastore
### Now you can create the Experimental Step describing your research experiment. Start by declaring an object of the ExpStep class, which you want to upload to the Datastore

#### An ExpStep Class instance looks like:

``` python
new_experiment = ExpStep(
    name: str = 'exp_name', # name of the experiment
    type: str = 'exp_type', # type of the experiment, has to be a defined sample type in the datastore
    metadata: dict = exp_metadata, # metadata of your data, is defined by the sample type
    space: str = 'exp_space', # space in which the experiment should be saved
    project: str = 'exp_project', # project in which the experiment should be saved
    collection: str = 'exp_space/exp_project/exp_collection', # collection in which the experiment should be saved
    parents: list[str] = [exp_parents], # parents of the experiment
    identifier: str = '', # identifier of the sample (SPACE/PROJECT/SAMPLE_CODE)
    permId: str = '', # permID of the sample
    sample_object=None, # Pybis object of the sample
    datasets: list = None, # dataset objects saved in the sample
    dataset_codes: list = None, # codes of the dataset objects saved in the sample
)
```

> Note: You do not have to define all of those attributes, only the name, type, metadata and the location of the Sample in the Database (space, project, collection) are necessary

### An example workflow which would upload your experimental step describing your concrete experiment could be 

#### 1. Define the ExpStep as you see fit
```python. 
from lebedigital.openbis.expstep import ExpStep

concrete_experiment = ExpStep()

concrete_experiment.name = 'Concrete_Experiment_1'
concrete_experiment.type = 'EXPERIMENTAL_STEP_CONCRETE_TEST'
concrete_experiment.metadata = {'name': 'Concrete_Experiment_1', 'date': '17:10:2022'}

# define the location in the datastore

concrete_experiment.space = 'SPACE'
concrete_experiment.project = 'PROJECT'
concrete_experiment.collection = 'SPACE/PROJECT/COLLECTION'
```

#### 2. Upload the ExpStep to the Database
```python
concrete_experiment.upload_expstep(o) # o = connected Interbis instance
```
#### Good to know:
* To upload a sample of a given type that type has to already exist in the datastore.
* The name of the sample has to be the same as the `$name` parameter in the metadata dictionary in order to be able to find the sample efficiently again.
* You can also load an existing sample into an ExpStep object using the `load_sample()` method of ExpStep.
* Every time you want to upload some Sample to the Database a Type Checker will run in order to make sure your ExpStep has been defined properly.
* You can look up which metadata the Experimental Step Type will accept by either getting the Excel import template running the Interbis `get_metadata_import_template()` method or get the properties of an Experimental Step Type by running the Interbis `get_sample_type_properties` method.

### You can now upload your datasets containing your raw or processed data files.

```python
concrete_experiment.upload_dataset(
    o, # A connected Interbis object
    props={
        '$name': 'dataset_name',
        'files': ['path_to_file'],
        'data_type': 'PROCESSED_DATA' # or RAW_DATA. PREVIEW, ... you can see them in the dataset types
    }
)
```

### If you have done everything correctly you should be able to see your Experimental Step and your Datasets in the Web-View of the Database in the location you defined in the ExpStep object. 