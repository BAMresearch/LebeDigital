import logging
from math import isnan
from sys import exit

import pandas as pd
import yaml
from termcolor import colored

from lebedigital.openbis.interbis import Interbis


class ExpStep:
    def __init__(
            self,
            name: str = '',
            type: str = '',
            metadata: dict = None,
            space: str = "",
            project: str = "",
            collection: str = "",
            parents: list = None,
            identifier: str = '',
            permId: str = '',
            sample_object=None,
            datasets: list = None,
            dataset_codes: list = None,
    ) -> None:
        """Creates an ExpStep Object as general python interface to the Experimental Step OpenBis object

        Args:
            name (str, optional): name of the experiment. Defaults to ''.
            type (str, optional): type of the experiment, has to be a defined sample type in the Datastore.
                Defaults to ''.
            metadata (dict, optional): metadata of your data, is defined by the sample type. Defaults to None.
            space (str, optional): space in which the experiment should be saved. Defaults to "".
            project (str, optional): project in which the experiment should be saved. Defaults to "".
            collection (str, optional): collection in which the experiment should be saved. Defaults to "".
            parents (list, optional): parents of the experiment. Defaults to None.
            identifier (str, optional): verbose identifier of the sample, indicates the space and project.
                Defaults to ''.
            permId (str, optional): PermID identifier of the sample. Defaults to ''.
            sample_object (Pybis, optional): the pybis sample object of the sample. Defaults to None.
            datasets (list, optional): list of dataset objects saved under the experimental step. Defaults to None.
            dataset_codes (list, optional): the PermIDs of the datasets. Defaults to None.
        """

        self.name = name
        self.type = type
        self.metadata = metadata
        self.space = space
        self.project = project
        self.collection = collection
        self.parents = parents
        self.identifier = identifier
        self.permId = permId
        self.sample_object = sample_object
        self.datasets = datasets
        self.dataset_codes = dataset_codes

        if not metadata:
            self.metadata = {}
        if not parents:
            self.parents = []
        if not datasets:
            self.datasets = []
        if not dataset_codes:
            self.dataset_codes = []

    def info(self):
        """Prints a description of the sample into the terminal/jupyter output
        """
        print('\n')
        title = colored(self.name, 'cyan')
        title = '=============== SAMPLE "' + title + '" ================'
        print(title)
        print('\n')
        for param, value in vars(self).items():

            column_param = colored(
                param, 'green') if value else colored(param, 'red')

            if isinstance(value, str):
                print(f'{column_param: >25}: {value}' + '\n')
            elif isinstance(value, dict):
                print(f'{column_param: >25}:')
                for key, val in value.items():
                    column_key = colored(key, 'green', 'on_grey')
                    print(' ' * 20 + f'{column_key}: {val}')
                print('\n')
            else:
                print(f'{column_param: >25}: {value}' + '\n')

        print('=' * int(float(len(title)) * 0.9))

    def import_from_template(self, path_to_file: str):
        """Imports the data from the template into the sample's metadata. To used with the generated import template
                function from InterBis which creates an Excel sheet for you to fill out.

        Args:
            path_to_file (str): Path to the filled out import template Excel sheet
        """

        df = pd.read_excel(path_to_file)
        df.columns = df.columns.str.lower()

        # if some params were left out then the row gets dropped
        df = df.dropna()

        # not useful for uploading to the datastore only for the user to know what to input
        df.pop('label')
        df.pop('description')

        # turn the df columns into a dict
        meta = dict(zip([par.lower() for par in df.param], df.value))
        self.metadata = meta

    # parse openbis into functions or a prepared dataframe, would need to update every upload script
    def check_type(self, o: Interbis):
        """Dry checking if the types of metadata entries are correct
           and if the upload step will run as expected

        Args:
            o (Interbis): currently running openbis interface

        Raises:
            ValueError: Type Mismatch between given and expected dtypes
            ValueError: NaN values in metadata
            ValueError: Unexpected parameters
            ValueError: Undeclared necessary attribute for upload
        """
        # TYPE CHECKING
        conv_dict = {
            'BOOLEAN': bool,
            'CONTROLLEDVOCABULARY': str,
            'DATE': str,
            'HYPERLINK': str,
            'INTEGER': int,
            'MATERIAL': str,
            'MULTILINE_VARCHAR': str,
            'OBJECT': str,
            'REAL': float,
            'TIMESTAMP': str,
            'VARCHAR': str,
            'XML': str,
            'SAMPLE': str
        }
        logging.debug('TYPE CHECKING OF THE SAMPLE')
        logging.debug('Setting up comparison dict')
        types_df = o.get_sample_type_properties(self.type)
        types_dict = dict(zip(types_df.code, types_df.dataType))
        logging.debug(types_dict)
        types_dict: dict = {k.lower(): conv_dict[v] for k, v in types_dict.items()}

        # Check if all keys are a subset of all possible keys
        logging.debug('Checking if keys are subset of defined parameters')
        if set(self.metadata.keys()).issubset(set(types_dict.keys())):
            for key, val in self.metadata.items():
                # Check if all values have the correct data type
                if not isinstance(val, types_dict[key]):
                    raise ValueError(
                        f'Type Checker: Type Mismatch.'
                        f'Entry "{key}" is of Type {type(val)} and should be {types_dict[key]}')
                    # Check if float values are not NaN
                if isinstance(val, float) and isnan(val):
                    raise ValueError(
                        f'Type Checker: NaN values are not accepted by Openbis.'
                        f'Entry "{key} has a value "{val}".')
        else:
            # Keys in metadata and not in object definition -> raises an error for every such key
            for key in self.metadata.keys():
                if key.lower() not in types_dict.keys():
                    raise ValueError(
                        f'Type Checker: Unexpected parameter.'
                        f'Key "{key}" is not in {self.type} defined metadata params')

        # Checking if all values are defined
        needed_attrs = ['name', 'type', 'collection', 'space', 'project']
        for attr, value in self.__dict__.items():
            if attr in needed_attrs:
                if not value:
                    raise ValueError(
                        f'Type Checker: Undeclared Attribute.'
                        f'Attribute "{attr}" has not been declared')

        # Printing warnings if parents are empty, does not break the function
        if not self.parents:
            logging.warning(
                'Type Checker: Parents attribute undeclared')

        # If you got to this line no game breaking bugs were found
        logging.info('Type Checker: Types are correctly set')

    def load_sample(self, o: Interbis, sample_identifier: str):
        """Fill the instance with properties of an experimental step in the datastore

        Args:
            o (Interbis): currently running openbis instance
            sample_identifier (str): identifier of the sample

        Returns:
            self: filled instance ExpStep
        """

        # Getting the sample properties from the Datastore
        sample_dict = o.get_sample_dict(sample_identifier)

        # Getting the list of properties only to filter them from sample_dict
        props_list = list(o.get_sample_type(sample_dict['type'])
                          .get_property_assignments()
                          .df['propertyType'])

        # Assigning the values from the dict to the new object
        self.name = sample_dict['$NAME']
        self.type = sample_dict['type']
        self.space = sample_dict['identifier'].split('/')[1]
        self.project = sample_dict['identifier'].split('/')[2]
        self.identifier = sample_dict['identifier']
        self.permId = sample_dict['permId']

        # Getting the object from the datastore to get the name of the collection and parents
        self.sample_object = o.get_sample(sample_dict['identifier'], props='*')
        self.collection = f'/{self.space}/{self.project}/{self.sample_object.experiment.code}'
        self.parents = self.sample_object.parents

        # Setting the metadata by filtering the sample_dict using props_list
        self.metadata = {key: sample_dict[key] for key in props_list}

        # Getting all datasets of the sample and their codes
        self.datasets = self.sample_object.get_datasets()
        self.dataset_codes = [ds.code for ds in self.datasets]

        return self

    def upload_expstep(self, o: Interbis, overwrite: bool = False) -> str:
        """Uploads the ExpStep object into the datastore

        Args:
            o (Interbis): Currently running openBIS interface
            overwrite (bool): Specifies if an existing sample should be overwritten. Defaults to False.

        Raises:
            ValueError: ValueError if the metadata is empty,
                        $name is mandatory.

        Returns:
            str: Name of the sample in the datastore
        """

        # TYPE CHECKING
        try:
            self.check_type(o)
        except Exception as e:
            logging.error(str(e))
            exit(1)

        # If a sample with the same name exists in the datastore you fetch it instead of creating a new one
        if o.exists_in_datastore(self.name):
            print(f'Sample {self.name} already exists in Datastore')
            samples_df = o.get_samples(
                type=self.type,
                space=self.space,
                project=self.project,
                props="$name"
            ).df

            # Getting the identifier from the dataframe
            samples_df.rename(columns={'$NAME': 'NAME'}, inplace=True)
            sample_identifier = samples_df.query(
                "NAME==@self.name")['identifier'].values[0]
            self.identifier = sample_identifier

            # Overwriting the sample by deleting and uploading a new one
            # TODO: Overwriting samples by changing their metadata instead of deleting and uploading new sample

            if overwrite:
                print('Overwriting the sample')
                self.delete_expstep(o, 'overwriting sample')
                sample = o.new_sample(
                    type=self.type,
                    space=self.space,
                    collection=self.collection,
                    parents=self.parents,
                    props=self.metadata,
                )
                sample.save()
                self.identifier = sample.identifier
            else:
                self.identifier = sample_identifier

        # No sample with the same name found -> uploading the sample
        else:
            print(f'Creating new sample {self.name}')
            sample = o.new_sample(
                type=self.type,
                space=self.space,
                collection=self.collection,
                parents=self.parents,
                props=self.metadata,
            )
            sample.save()
            self.identifier = sample.identifier
            self.identifier = sample.identifier
        return self.identifier

    def delete_expstep(self, o: Interbis, reason: str):
        """Deletes a sample from the datastore

        Args:
            o (Interbis): Currently running openbis interface
            reason (str): Reason for deletion
        """
        # Check if the sample exists and delete it
        if o.exists_in_datastore(self.name):
            o.get_sample(self.identifier).delete(reason)
        else:
            raise ValueError('Sample not in datastore')

    def upload_dataset(self, o: Interbis, props: dict):
        """Uploads a dataset to the ExpStep.
        
        Requires a dictionary with $name, files and data_type

        Args:
            o (Interbis): Currently running openbis interface
            props (str, optional): Metadata of the dataset.

        Returns:
            str: Properties of the dataset
        """

        # Checking if the name of the dataset is included in props
        if '$name' in props:
            exists_test = o.get_datasets(where={'$name': props['$name']})
        elif '$NAME' in props:
            exists_test = o.get_datasets(where={'$NAME': props['$NAME']})
        else:
            raise KeyError('$name not defined in props')

        # Checking if the files and data_type are specified in props
        if 'files' not in props:
            raise KeyError(f'files not specified in props')
        if 'data_type' not in props:
            raise KeyError(f'data_type not specified in props')

        files = props.pop('files')
        data_type = props.pop('data_type')

        # If a dataset with the same name was found in the datastore that dataset will be returned and none will be
        # uploaded
        if exists_test:
            name = props['$name']
            print(
                f'Dataset(s) with the same name already present in the Datastore.\nTo upload the dataset you must '
                f'first delete the other dataset with name {name}')
            ds = o.get_datasets(where={'$name': props['$name']})

        # Uploading the dataset
        else:
            ds = o.new_dataset(
                type=data_type,
                collection=self.collection,
                sample=self.identifier,
                files=files,
                props=props
            )
            ds.save()

        return ds

    def download_datasets(self, o: Interbis, path: str, data_type: str = ''):
        """Downloads all datasets which are assigned to that sample

        Args:
            o (Interbis): Currently running openBIS interface
            path (str): Path where the datasets should be saved
            data_type (str): If specified will only download data sets of that type

        Raises:
            ValueError: Raises an error when no datasets are found under the sample
        """
        if not o.is_session_active():
            raise ValueError('Interbis object is not connected')

        if not data_type:
            length = len(self.datasets)
        else:
            # We get all data_types from the dataset object list and filter them using the given data_type
            # the resulting list will have the correct length
            ds_list = [s.type.code for s in self.datasets]
            ds_list_filtered = list(filter(
                lambda x: x == data_type, ds_list
            ))
            length = len(ds_list_filtered)

        # If the sample has no datasets (of given data_type) an error will be thrown
        if not length:
            err_msg = "No Datasets found under the sample" if not data_type else f"No Datasets of type {data_type} " \
                                                                                 f"found under the sample"
            raise ValueError(err_msg)

        file_plural = 'FILES' if length > 1 else 'FILE-'

        print(
            f'----------DOWNLOADING {length} {file_plural}----------\n')
        for dataset in self.datasets:
            # If data_type was specified download only the datasets with that data type
            if data_type:
                if dataset.type == data_type:
                    print(f'Downloading dataset {dataset.code}')
                    print(f'Files: {dataset.file_list}')
                    dataset.download(
                        destination=path,
                        create_default_folders=False,
                        wait_until_finished=False,
                    )
                    print("\n")

            # If data_type was NOT specified download all datasets
            else:
                print(f'Downloading dataset {dataset.code}')
                print(f'Files: {dataset.file_list}\n')
                dataset.download(
                    destination=path,
                    create_default_folders=False,
                    wait_until_finished=False,
                )
                print("\n")

        print('----------DOWNLOAD FINISHED----------')

    def save_expstep_yaml(self, yaml_path: str):
        """Saves a log file of a sample in openBIS. Used in doit upload

        Args:
            yaml_path (str): Path where the yaml file should be saved
        """

        modified_dict = self.__dict__

        # The python objects are ignored in the printout, too long and too useless
        modified_dict['sample_object'] = '--A PYTHON OBJECT WAS HERE--'
        modified_dict['datasets'] = '--A PYTHON OBJECT WAS HERE--'

        with open(yaml_path, 'w') as file:
            _ = yaml.dump(modified_dict, file)

    def save_sample_yaml(self, o: Interbis, yaml_path: str):
        df = o.get_sample_type_properties(self.type)

        # deleting unnecessary columns
        df = df.drop(df.iloc[:, 5:], axis=1)
        df = df.drop(df.columns[3], axis=1)

        # reordering the columns to fit the predefined order
        df = df.iloc[:, [0, 3, 1, 2]]

        # setting the index to the future dict key
        df = df.set_index('code')

        # making the output dict, filtering the keys from the output dict
        out = df.to_dict("index")
        out = {key: [*val.values()] for key, val in out.items()}

        with open(yaml_path, 'w') as file:
            _ = yaml.dump(out, file)
