# Interbis/ExpStep

## Pybis Extension for improving data synchronisation with openBIS

## Example use case of the ExpStep class

### The main part of the package is the ExpStep class. Start by declaring an object of the ExpStep class, which you want to upload to the Datastore

``` python
new_experiment = ExpStep(
    name: str = exp_name, # name of the experiment
    type: str = exp_type, # type of the experiment, has to be a defined sample type in the Datastore
    metadata: dict = exp_metadata, # metadata of your data, is defined by the sample type
    space: str = exp_space, # space in which the experiment should be saved
    project: str = exp_project, # project in which the experiment should be saved
    collection: str = exp_collection, # collection in which the expeiment should be saved
    parents: str = exp_parents, # parents of the experiment
    identifier: str = '', # identifier of te sample (SPACE/PROJECT/SAMPLE_CODE)
    permId: str = '', # permID of the sample
    sample_object=None, # Pybis object of the sample
    datasets: list = None, # Dataset objects saved in the Sample
    dataset_codes: list = None, # Codes of the Dataset objects saved in the sample
)
```

### In order to use ExpStep, first create an expstep instance of the class, then set the attributes of your sample like name, type and metadata
```python
from lebedigital.openbis.expstep import ExpStep

my_experiment = ExpStep()

my_experiment.name = 'Experiment 1'
my_experiment.type = 'EXPERIMENTAL_STEP_MYTYPE'
my_experiment.metadata = {'param1': 'val1', 'param2': 'val2'}
...
```

> Then you can log into your datastore and upload the sample

```python
from lebedigital.openbis.interbis import Interbis

o = Interbis('https://yoururl.com')
o.connect_to_datastore()

my_experiment.upload_expstep(o)
```