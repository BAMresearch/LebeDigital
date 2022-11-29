# Interbis/ExpStep

## Pybis Extension for improving data synchronisation with openBIS

## Example use case of the ExpStep class

### The main part of the package is the ExpStep class. Start by declaring an object of the ExpStep class, which you want to upload to the Datastore

``` python
new_experiment = ExpStep(
    name: str = exp_name, # name of the experiment
    type: str = exp_type, # type of the experiment, has to be a defined sample type in the datastore
    metadata: dict = exp_metadata, # metadata of your data, is defined by the sample type
    space: str = exp_space, # space in which the experiment should be saved
    project: str = exp_project, # project in which the experiment should be saved
    collection: str = exp_collection, # collection in which the experiment should be saved
    parents: str = exp_parents, # parents of the experiment
    identifier: str = '', # identifier of the sample (SPACE/PROJECT/SAMPLE_CODE)
    permId: str = '', # permID of the sample
    sample_object=None, # Pybis object of the sample
    datasets: list = None, # dataset objects saved in the sample
    dataset_codes: list = None, # codes of the dataset objects saved in the sample
)
```

### In order to use ExpStep, first create an ExpStep instance of the class, then set the attributes of your sample like name, type and metadata
* You also have to specify the location in the datastore
```python. 
from lebedigital.openbis.expstep import ExpStep

my_experiment = ExpStep()

my_experiment.name = 'Experiment 1'
my_experiment.type = 'EXPERIMENTAL_STEP_MYTYPE'
my_experiment.metadata = {'param1': 'val1', 'param2': 'val2'}

# location in the datastore

my_experiment.space = 'SPACE'
my_experiment.project = 'PROJECT'
my_experiment.collection = 'SPACE/PROJECT/COLLECTION'
...
```

* Then you can log into your datastore and upload the sample using the openbis interface class Interbis

```python
from lebedigital.openbis.interbis import Interbis

o = Interbis('https://yoururl.com')
o.connect_to_datastore()

my_experiment.upload_expstep(o)
```
* To upload a sample of a given type that type has to already exist in the datastore. You can define custom types by using the ELN Web View or do it from the notebook/script by using the o.create_sample_type() function.

```python
o.create_sample_type(
    sample_code='EXPERIMENTAL_STEP_MYTYPE',
    sample_prefix='MYTYPE',
    sample_properties={
        'param1': ['VARCHAR', 'my_label1', 'my_description1'],
        'param2': ['VARCHAR', 'my_label2', 'my_description2'],
    }
)
```
* You can upload datasets by using the upload_dataset() function.

```python
my_sample.upload_dataset(
    o, # A connected Interbis object
    props={
        '$name': 'dataset_name',
        'files': ['path_to_file'],
        'data_type': 'PROCESSED_DATA' # or RAW_DATA. PREVIEW, ... you can see them in the dataset types
    }
)
```