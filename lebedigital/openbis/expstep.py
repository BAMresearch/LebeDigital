from math import isnan
from pybis import Openbis
import pandas as pd
import json
import os
import logging
from getpass import getpass
import sys
from sys import exit


class ExpStep:
    def __init__(
            self,
            name: str = '',
            type: str = '',
            data_path: list = [],
            data_type: str = '',
            code: str = '',
            metadata: dict = {},
            space: str = "",
            project: str = "",
            collection: str = "",
            parents: list = [],
            children: list = [],
            identifier: str = '',
            permId: str = '',
            sample_object=None,
            datasets: list = [],
            dataset_codes: list = [],
    ) -> None:
        """Creates an ExpStep Object

        Args:
            name (str): name of the experiment
            type (str): type of the experiment, has to be a defined sample type in the Datastore
            data_path (str, optional): local path to the data file. Defaults to ''.
            data_type (str, optional): type of the data, has to be a defined data type in the Datastore. Defaults to ''.
            code (str, optional): the code the experiment will have in the Datastore. Defaults to ''.
            metadata (dict, optional): metadata of your data, is defined by the sample type. Defaults to {}.
            space (str, optional): space in which the experiment should be saved. Defaults to "".
            project (str, optional): project in which the experiment should be saved. Defaults to "".
            collection (str, optional): collection in which the expeiment should be saved. Defaults to "".
            parents (list, optional): parents of the experiment. Defaults to [].
            children (list, optional): children of the experiment. Defaults to [].
        """

        self.name = name
        self.type = type
        self.data_path = data_path
        self.data_type = data_type
        self.metadata = metadata
        self.space = space
        self.project = project
        self.collection = collection
        self.parents = parents
        self.children = children
        self.code = code
        self.identifier = identifier
        self.permId = permId
        self.sample_object = sample_object
        self.datasets = datasets
        self.dataset_codes = dataset_codes

    @staticmethod
    def connect_to_datastore(url='https://test.datastore.bam.de/openbis/', *args, **kwargs):
        """Connects to a Datastore

        Args:
            url (str, optional): Datastore URL.
                    Defaults to 'https://test.datastore.bam.de/openbis/'.

        Returns:
            Openbis: Connected and running Datastore instance
        """
        os.environ['OPENBIS_HOST'] = url
        o = Openbis(os.environ['OPENBIS_HOST'])

        # If a session is already active (when running a notebook again) return the active object
        if o.is_session_active():
            return o

        # We get the username from the os username and the password has to be given by the user at runtime
        else:
            os.environ['OPENBIS_USERNAME'] = os.getlogin()
            if 'password' in kwargs:
                os.environ['OPENBIS_PASSWORD'] = kwargs['password']
            else:
                os.environ['OPENBIS_PASSWORD'] = getpass("Give Password: ")
            try:
                o.login(os.environ['OPENBIS_USERNAME'],
                        os.environ['OPENBIS_PASSWORD'])
            except ValueError:
                print("Wrong Credentials")
                exit(1)

            return o

    @staticmethod
    def read_metadata_json(json_path: str) -> dict:
        """Read the metadata from a json file. You can download the metadata as json from the datastore website.

        Args:
            json_path (str): Path to json file

        Returns:
            dict: Json file as dict. Can be directly imported into the
                  datastore.
        """
        with open(json_path, 'r') as f:
            metadata = json.load(f).get('properties')
        # making keys lowercase, k - key, v - value
        metadata = {k.lower(): v for k, v in metadata.items()}
        return metadata

    @staticmethod
    def gen_metadata_import_template(
        o: Openbis,
        sample_type: str,
        write: bool = False,
        sheet_name: str = 'metadata',
        path: str = '',
    ):
        """ Generate an excel data import template for a given object type.

        Args:
            o (Openbis): Currently running openbis instance
            sample_type (str): The type of the sample (ex. EXPERIMENTAL_STEP_SOMETHING)
            write (bool, optional): Specifies if the function should return the template as a pandas DataFrame or if it should be written to an excel file. Defaults to False.
            sheet_name (str, optional): The name of the sheet in an excel file. Defaults to 'metadata'.
            path (str, optional): Path where the excel sheet should be saved. Only use with write set to True. Defaults to ''.

        Raises:
            ValueError: Raises an error when write is set to True and no path is gven

        Returns:
            pandas.DataFrame: returns a DataFrame with the import template 
        """

        # Raises an error when no path given and write set to True
        if write and not path:
            raise ValueError("No path given")

        # Get the properties of the sample type from the datastore
        # and turn it into a dict

        meta_list = list(o.get_sample_type(
            sample_type).get_property_assignments().df['propertyType'])
        meta_dict = {idx: v for idx, v in enumerate(meta_list)}

        # Parse the dict into a dataframe

        meta_df = pd.DataFrame.from_dict(
            meta_dict, orient="index", columns=['Param'])

        # Get property descriptions
        props_list = list(o.get_sample_type(sample_type)
                          .get_property_assignments()
                          .df['propertyType'])

        # Build a dataframe with the properties and their descriptions
        df = pd.DataFrame()
        for prop in props_list:
            pt = o.get_property_type(prop)
            prop_series = pd.Series(pt.attrs.all())
            prop_df = pd.DataFrame(
                {prop: prop_series.values}, index=prop_series.index)

            if df.empty:
                df = prop_df
            else:
                df = df.join(prop_df)

        prop_df = df

        meta_df['Label'] = prop_df.loc['label'].values
        meta_df['Description'] = prop_df.loc['description'].values

        # Add an empty column where you can input the values
        meta_df = meta_df.assign(Value='')

        # If write was specified then write to file
        if write:
            # Check whether sheet already exists, if yes set args as to replace
            # the current sheet with the new one
            (mode, args) = ('a', {'if_sheet_exists': 'replace'}) if os.path.exists(
                path) else ('w', {})

            # Write the sheet
            with pd.ExcelWriter(path=path, mode=mode, **args) as writer:
                meta_df.to_excel(excel_writer=writer, sheet_name=sheet_name)

        else:
            return meta_df

    # The simple static methods will get integrated into their methods later because they are useless to define over pybis itself
    @staticmethod
    def get_space_names(o: Openbis):
        spaces_df = o.get_spaces().df
        return list(spaces_df['code'].values)

    @staticmethod
    def get_project_names(o: Openbis, space: str):
        projects_df = o.get_projects(space=space).df
        return [name.split('/')[-1] for name in list(projects_df['identifier'].values)]

    @staticmethod
    def get_experiment_names(o: Openbis, space: str, project: str):
        experiments_df = o.get_experiments(
            space=space,
            project=project,
        ).df
        return [name.split('/')[-1] for name in list(experiments_df['identifier'].values)]

    get_collection_names = get_experiment_names

    @staticmethod
    def get_samples(o: Openbis, space: str, project: str, collection: str, **kwargs):

        if "experiment" in kwargs:
            kwargs["collection"] = kwargs["experiment"]
            kwargs.pop("experiment", None)

        return o.get_samples(space=space, project=project, collection=collection)

    get_objects = get_samples

    @staticmethod
    def get_sample_names(o: Openbis, space: str, project: str, collection: str, **kwargs):

        if "experiment" in kwargs:
            kwargs["collection"] = kwargs["experiment"]
            kwargs.pop("experiment", None)

        samples_df = o.get_samples(
            space=space,
            project=project,
            collection=collection,
            props=['$name']
        ).df
        samples_df.columns = samples_df.columns.str.upper()
        sample_codes = [name.split('/')[-1]
                        for name in list(samples_df['IDENTIFIER'].values)]
        sample_names = list(samples_df['$NAME'].values)
        return [f'{code} ({name})' for code, name in zip(sample_codes, sample_names)]

    get_object_names = get_sample_names

    @staticmethod
    def get_sample_codes(o: Openbis, space: str, project: str, collection: str, **kwargs):

        if "experiment" in kwargs:
            kwargs["collection"] = kwargs["experiment"]
            kwargs.pop("experiment", None)

        samples_df = o.get_samples(
            space=space,
            project=project,
            collection=collection,
            props=['$name']
        ).df
        return [name.split('/')[-1] for name in list(samples_df['identifier'].values)]

    @staticmethod
    def get_sample_dict(o: Openbis, identifier) -> dict:
        """ Returns all/some metadata about the sample

        Args:
            o (Openbis): Currently running openbis instance
            identifier (_type_): Identifier of the sample

        Returns:
            dict: Dict containing the information as metadata_label : metadata_value
        """
        sample = o.get_sample(identifier)
        sample_info_dict = sample.attrs.all()
        sample_prop_dict = {key.upper(): val for key,
                            val in sample.p.all().items()}
        return sample_info_dict | sample_prop_dict

    @staticmethod
    def get_overview(o: Openbis, level: str, **kwargs) -> dict:
        """ Generates an overview for the samples stored in the datastore
            You need to provide the 
        Args:
            o (Openbis): Currently running openbis instance
            level (str): What entity should be the highest level of the overview (space/project/collection overview)

        Raises:
            ValueError: Raises an error when no correct value for the level was specified

        Returns:
            dict: Returns a dictionary with the overview
        """

        # We check if the openbis aliases were specified in the function arguments
        if "experiment" in kwargs:
            kwargs["collection"] = kwargs["experiment"]
            kwargs.pop("experiment", None)

        if "object" in kwargs:
            kwargs["sample"] = kwargs["object"]
            kwargs.pop("object", None)

        # We go through all entries in for loops, the only difference between levels is where we start the loop
        if level == 'full':
            space_dict = {'DATASTORE': {}}

            for space in ExpStep.get_space_names(o):
                project_dict = {}
                for project in ExpStep.get_project_names(o, space):
                    collection_dict = {}
                    for collection in ExpStep.get_collection_names(o, space, project):
                        sample_list = []
                        for sample in ExpStep.get_sample_names(o, space, project, collection):
                            sample_list.append(sample)
                        collection_dict[collection] = sample_list
                    project_dict[project] = collection_dict
                space_dict['DATASTORE'][space] = project_dict

            return space_dict

        elif level == 'space':
            space = kwargs.pop('space')

            # project_dict = {space: {}}
            project_dict = {space: {}}
            for project in ExpStep.get_project_names(o, space):
                collection_dict = {}
                for collection in ExpStep.get_collection_names(o, space, project):
                    sample_list = []
                    for sample in ExpStep.get_sample_names(o, space, project, collection):
                        sample_list.append(sample)
                    collection_dict[collection] = sample_list

                project_dict[space][project] = collection_dict

            return project_dict

        elif level == 'project':
            space = kwargs.pop('space')
            project = kwargs.pop('project')

            collection_dict = {project: {}}
            for collection in ExpStep.get_collection_names(o, space, project):
                sample_list = []
                for sample in ExpStep.get_sample_names(o, space, project, collection):
                    sample_list.append(sample)
                collection_dict[project][collection] = sample_list

            return collection_dict

        elif level == 'collection':
            space = kwargs.pop('space')
            project = kwargs.pop('project')
            collection = kwargs.pop('collection')

            sample_list = []
            for sample in ExpStep.get_sample_names(o, space, project, collection):
                sample_list.append(sample)
            sample_dict = {collection: sample_list}

            return sample_dict

        else:
            raise ValueError('No correct level specified')

    @staticmethod
    def get_props_list(o: Openbis, identifier: str) -> list:
        return [prop.upper() for prop in list(o.get_sample(identifier).p().keys())]

    def sync_name(self):
        """ Synchronizes the name of the attribute with the name in metadata.
            Have to be the same in order to search for the sample accurately
        """
        new_name = getattr(self, "metadata").get('$name')
        setattr(self, "name", new_name)

    def sync_code(self):
        """ Synchronizes code with name using the correct code format.
            Deprecated, dont use will be removed lated.
        """
        name = self.name
        allowed_chars = "1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-_"
        name = "".join(list(filter(lambda x: x != ',', name)))

        code = ''
        for char_idx in range(len(name)):
            if not name[char_idx] in allowed_chars:
                code += '_'
            else:
                code += name[char_idx]

        setattr(self, 'code', code)

    def metadata_import_template(
            self,
            o: Openbis,
            path: str = '',
            sheet_name: str = 'metadata',
            write: bool = False):
        """ Generates an excel sheet to be filled out by hand and then
            imported into an experimental step

        Args:
            o (Openbis): Currently running openbis instance
            path (str): Path where the sheet should be saved.
                        Must include the file name and extenstion
            sheet_name (str, optional): Name of the sheet.
                                        Defaults to 'metadata'.
        """
        if write and not path:
            raise ValueError("No path given")

        # Get the properties of the sample type from the datastore
        # and turn it into a dict

        meta_list = list(o.get_sample_type(
            self.type).get_property_assignments().df['propertyType'])
        meta_dict = {idx: v for idx, v in enumerate(meta_list)}

        # Parse the dict into a dataframe
        meta_df = pd.DataFrame.from_dict(
            meta_dict, orient="index", columns=['Param'])

        # Get property descriptions
        prop_df = self.get_property_types(o)
        meta_df['Description'] = prop_df.loc['label'].values

        # Add an empty column where you can input the values
        meta_df = meta_df.assign(Value='')

        # If write was specified then write to file
        if write:
            # Check whether sheet already exists, if yes set args as to replace
            # the current sheet with the new one
            (mode, args) = ('a', {'if_sheet_exists': 'replace'}) if os.path.exists(
                path) else ('w', {})

            # Write the sheet
            with pd.ExcelWriter(path=path, mode=mode, **args) as writer:
                meta_df.to_excel(excel_writer=writer, sheet_name=sheet_name)

        else:
            return meta_df

    def exists_in_datastore(self, o: Openbis) -> bool:
        """Checks whether a sample is already in the datastore

        Args:
            openbis (Openbis): currently running openbis instance

        Returns:
            bool: Outcome of the check
        """
        samples = o.get_samples(
            space=self.space, project=self.project, props="$name")
        samples_df = samples.df
        samples_df.columns = samples_df.columns.str.upper()
        return self.name in samples_df['$NAME'].values

    def find_collection(
        self,
        o: Openbis,
        collection_name: str,
        id_type: int = 0,
    ) -> str:
        """Finds a collection in the Datastore

        Args:
            o (Openbis): Currently running Openbis instance
            collection_name (str): Name of the Collection
            id_type (int, optional): Specifies return type to identifier or permId. Defaults to 0 = permId.

        Returns:
            str: identifier or permId of the Collection
        """

        id_val = 'identifier' if id_type else 'permId'

        collections_df = o.get_space(self.space).get_project(
            self.project).get_experiments().df
        return collections_df[collections_df['identifier'].str.contains(
            collection_name)].iloc[0][id_val]

    def check_type(self, o: Openbis):
        """Dry checking if the types of metadata entries are correct
           and if the upload step will run as expected

        Args:
            o (Openbis): currently running openbis instance

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

        types_df = self.get_property_types(o)
        types_dict = types_df['dataType'].to_dict()
        types_dict = {k.lower(): conv_dict[v] for k, v in types_dict.items()}

        # Check if all keys are a subset of all possible keys
        if set(self.metadata.keys()).issubset(set(types_dict.keys())):
            for key, val in self.metadata.items():
                # Check if all values have the correct data type
                if not isinstance(val, types_dict[key]):
                    raise ValueError(
                        f'Type Mismatch. Entry "{key}" is of Type {type(val)} and should be {types_dict[key]}')
                # Check if float values are not NaN
                if isinstance(val, float) and isnan(val):
                    raise ValueError(
                        f'NaN values are not accepted by Openbis. Entry "{key} has a value "{val}".')
        else:
            # Keys in metadata and not in object definition -> raises an error for every such key
            for key in self.metadata.keys():
                if key not in types_dict.keys():
                    raise ValueError(
                        f'Unexpected parameter. Key "{key}" is not in {self.type} defined metadata params')

        # Checking if all values are defined
        needed_attrs = ['name', 'type', 'collection', 'space', 'project']
        for attr, value in self.__dict__.items():
            if attr in needed_attrs:
                if not value:
                    raise ValueError(
                        f'Undeclared Attribute. Attribute "{attr}" has not been delcared')

        # Printing warnings if parents are empty, does not break the function
        if not self.parents:
            logging.warning(
                'parents attribute undeclared')

        # If you got to this line no gamebreaking bugs were found
        print('Types are correctly set')

    def upload_expstep(self, o: Openbis, overwrite: bool = False) -> str:
        """Uploads the ExpStep object into the datastore

        Args:
            o (Openbis): Currently running openbis instance
            name (str): Name of the sample
            overwrite (bool): Specifies if an existing Sample should be overwritten. Defaults to False.

        Raises:
            ValueError: ValueError if the Metadata is empty,
                        $name is mandatory.

        Returns:
            str: Name of the sample in the datastore
        """

        self.sync_name()
        # self.sync_code()

        # TYPE CHECKING
        try:
            self.check_type(o)
        except Exception as e:
            logging.error(str(e))
            exit(1)

        # If a sample with the same name exists in the Datastore you fetch it instead of creating a new one
        if self.exists_in_datastore(o):
            print(f'Sample {self.name} already exists in Datastore')
            samples_df = o.get_samples(
                type=self.type,
                space=self.space,
                project=self.project,
                props="$name"
            ).df

            # Gettingthe identifier from the dataframe
            # TODO: Use pybis queries (where) here instead of pulling the entire dataframe

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
                    # code=self.code,
                    type=self.type,
                    space=self.space,
                    collection=self.collection,
                    parents=self.parents,
                    children=self.children,
                    props=self.metadata,
                )
                sample.save()
                return sample.identifier
            else:
                return sample_identifier

        # No sample with the same name found -> uploading the sample
        else:
            print(f'Creating new sample {self.name}')
            sample = o.new_sample(
                # code=self.code,
                type=self.type,
                space=self.space,
                collection=self.collection,
                parents=self.parents,
                children=self.children,
                props=self.metadata,
            )
            sample.save()
            self.identifier = sample.identifier
            return sample.identifier

    @staticmethod
    def load_sample(o: Openbis, sample_identifier: str):
        """Loads an expstep from an experimental step in the datastore with its properties

        Args:
            o (Openbis): currently running openbis instance
            sample_identifier (dict): identifier of the sample

        Returns:
            ExpStep: ExpStep object containing metadata of the sample
        """

        # Getting the properties of the sample
        sample_dict = ExpStep.get_sample_dict(o, sample_identifier)

        #  Getting a list of the properties only to flter them from the sample_dict
        props_list = list(o.get_sample_type(sample_dict['type'])
                          .get_property_assignments()
                          .df['propertyType'])

        # Getting the name of the sample
        sample_name = sample_dict['identifier'].split('/')[3]

        # Getting the sample
        sample = o.get_sample(sample_dict['identifier'], props='*')

        # Getting the name of the collection
        sample_collection = sample.experiment.code

        # Getting the name of the parents
        sample_parents = sample.parents

        # Getting the metadata from the sample by comparing it with the list of the properties
        sample_metadata = dict((key, sample_dict[key]) for key in props_list)

        # Getting the datasets uploaded to the sample
        sample_datasets = sample.get_datasets()

        # Getting the dataset codes
        sample_dataset_codes = [ds.code for ds in sample_datasets]

        # Combining all together to build an ExpStep object
        sample_step = ExpStep(
            name=sample_dict['$NAME'],
            type=sample_dict['type'],
            space=sample_dict['identifier'].split('/')[1],
            project=sample_dict['identifier'].split('/')[2],
            collection=sample_collection,
            parents=sample_parents,
            metadata=sample_metadata,
            code=sample_name,
            identifier=sample_dict['identifier'],
            permId=sample_dict['permId'],
            sample_object=sample,
            datasets=sample_datasets,
            dataset_codes=sample_dataset_codes
        )

        return sample_step

    def fill_sample(self, o: Openbis, sample_identifier: str):
        """Fill the instance with properties of an experimental step in the datastore

        Args:
            o (Openbis): currently running openbis instance
            sample_identifier (str): identifier of the sample

        Returns:
            self: filled instance ExpStep
        """

        # Getting the sample properties from the Datastore
        sample_dict = self.get_sample_dict(o, sample_identifier)

        # Getting the list of properties only to filter them from sample_dict
        props_list = list(o.get_sample_type(sample_dict['type'])
                          .get_property_assignments()
                          .df['propertyType'])

        # Assigning the values from the dict to the new object
        self.name = sample_dict['$NAME']
        self.type = sample_dict['type'],
        self.space, = sample_dict['identifier'].split('/')[1],
        self.project, = sample_dict['identifier'].split('/')[2],
        self.identifier = sample_dict['identifier']
        self.permId = sample_dict['permId']
        self.code = sample_dict['identifier'].split('/')[3]

        # Getting the object from the datastore to get the name of the collection and parents
        self.sample_object = o.get_sample(sample_dict['identifier'], props='*')
        self.collection = f'/{self.space}/{self.project}/{self.sample_object.experiment.code}'
        self.parents = self.sample_object.parents

        # Setting the metadata by filtering the sample_dict using props_list
        self.metadata = {key: sample_dict[key] for key in props_list}

        # Getting all datasets of the sample and their codes
        self.datasets = self.sample_object.get_datasets()
        self.dataset_codes = [ds.code for ds in self.datasets]

        self.data_path = []
        self.data_type = ''
        self.children = []

        return self

    def get_property_types(self, o: Openbis) -> pd.DataFrame:
        """Returns a DataFrame of the sample properties with their descriptions, labels and other metadata

        Args:
            o (Openbis): currently running openbis datastore

        Returns:
            pd.DataFrame: DataFrame of all properties with their attributes
        """

        # Getting a list of all the samples properties
        props_list = list(o.get_sample_type(self.type)
                          .get_property_assignments()
                          .df['propertyType'])

        df = pd.DataFrame()
        # Getting the metadata of every entry in props_list
        for prop in props_list:
            pt = o.get_property_type(prop)
            prop_series = pd.Series(pt.attrs.all())
            prop_df = pd.DataFrame(
                {prop: prop_series.values}, index=prop_series.index)

            # Combining the props together into a dataframe
            if df.empty:
                df = prop_df
            else:
                df = df.join(prop_df)

        return df.transpose()

    def delete_expstep(self, o: Openbis, reason: str):
        """Deletes a sample from the datastore

        Args:
            o (Openbis): Currently running openbis datastore
            reason (str): Reason for deletion
        """
        # Check if the sample exsts and delete it
        if self.exists_in_datastore(o):
            o.get_sample(self.identifier).delete(reason)
        else:
            raise ValueError('Sample not in datastore')

    def upload_dataset(self, o: Openbis, props: dict):
        """Uploads a dataset to the ExpStep

        Args:
            o (Openbis): Currently running openbis datastore
            props (str, optional): Metadata of the dataset.

        Returns:
            str: Properties of the dataset
        """

        # Checking if the name of the dataset is included in props
        test = ''
        if '$name' in props:
            test = o.get_datasets(where={'$name': props['$name']})
        elif '$NAME' in props:
            test = o.get_datasets(where={'$NAME': props['$NAME']})
        else:
            raise ValueError('$name not defined in props')

        # Checking if the files and data_type are specified in props
        if not 'files' in props:
            raise KeyError(f'files not specified in props')
        if not 'data_type' in props:
            raise KeyError(f'data_type not specified in props')

        files = props.pop('files')
        data_type = props.pop('data_type')

        # If a dataset with the same name was found in the datastore that dataset will be returned and none will be uploaded
        if test:
            name = props['$name']
            print(
                f'Dataset(s) with the same name already present in the Datastore.\nTo upload the dataset you must first delete the other dataset with name {name}')
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

    def download_datasets(self, o: Openbis, path: str):
        """Downloads all datasets which are asigned to that sample
            TODO: Add option to specify which dataset type should be downloaded

        Args:
            o (Openbis): Currently running openbis instance
            path (str): path where the datasets should be saved

        Raises:
            ValueError: Raises an error when no datasets are found under the sample
        """

        # If the sample has no datasets an error will be thrown
        if not len(self.datasets):
            raise ValueError('No Datasets found under the sample')

        file_plural = 'FILES' if len(self.datasets) > 1 else 'FILE-'

        self.data_path = []

        print(
            f'----------DOWNLOADING {len(self.datasets)} {file_plural}----------')
        for dataset in self.datasets:
            print(f'Downloading dataset {dataset.code}')
            print(f'Files: {dataset.file_list}')
            dataset.download(
                destination=path,
                create_default_folders=False,
                wait_until_finished=False,
            )
            if self.data_type:
                print(
                    'Warning: Data type was already specified before, may cause issues when uploading data')
                print(f'Old data type: {self.data_type}')
                self.data_type = dataset.type.code
                print(f'New data type: {self.data_type}')
            else:
                self.data_type = dataset.type.code

            self.data_path.append(f'{path}/{dataset}')
        print('----------DOWNLOADING FINISHED----------')

    def create_sample_type(self, o: Openbis, sample_code: str, sample_prefix: str, sample_properties: list[str]):
        """Creating a new sample type to be used by samples/objects in the datastore

        Args:
            o (Openbis): A currently running Openbis instance
        """

        # SUPRESSING PRINTS FOR ASSIGNING PROPERTIES
        # Disable
        def blockPrint():
            sys.stdout = open(os.devnull, 'w')

        # Restore
        def enablePrint():
            sys.stdout = sys.__stdout__

        sample_types = o.get_sample_types().df
        if sample_code in list(sample_types['code']):
            print(f'Sample type {sample_code} already exists in the Datastore')

        else:
            new_sample_type = o.new_sample_type(
                code=sample_code,
                generatedCodePrefix=sample_prefix,
                # description='Testing Experimental Step', DOES NOT WORK WITH PYBIS, PYBIS DOES NOT ACCEPT DESCTIPTION ARGUMENT
                autoGeneratedCode=True,
                subcodeUnique=False,
                listable=True,
                showContainer=False,
                showParents=True,
                showParentMetadata=True,
                validationPlugin='EXPERIMENTAL_STEP.date_range_validation'
            )
            new_sample_type.save()
            print(f'Sample type {sample_code} created')

        return o.get_sample_type(sample_code)


def new_object_text():

    # SUPRESSING PRINTS FOR ASSIGNING PROPERTIES
    # Disable
    def blockPrint():
        sys.stdout = open(os.devnull, 'w')

    # Restore
    def enablePrint():
        sys.stdout = sys.__stdout__

    o = ExpStep.connect_to_datastore()

    sample_code = 'EXPERIMENTAL_STEP_TEST'
    sample_types = o.get_sample_types().df
    created_sample = False

    # Creating new sample type or fetching an existing one if it exists

    if not sample_code in list(sample_types['code']):

        test_sample_type = o.new_sample_type(
            code=sample_code,
            generatedCodePrefix='TEST',
            # description         = 'Testing Experimental Step', DOES NOT WORK WITH PYBIS, PYBIS DOES NOT ACCEPT DESCTIPTION ARGUMENT
            autoGeneratedCode=True,
            subcodeUnique=False,
            listable=True,
            showContainer=False,
            showParents=True,
            showParentMetadata=True,
            validationPlugin='EXPERIMENTAL_STEP.date_range_validation'
        )
        test_sample_type.save()
        created_sample = True
        print(f'Creating new sample type {test_sample_type}')

    else:
        test_sample_type = o.get_sample_type(sample_code)
        print(f'Fetching existing sample type {test_sample_type}')

    # CREATING NEW VOCABULARY

    voc_list = list(o.get_vocabularies().df['code'])
    created_voc = False
    voc_code = 'TEST_VOC'

    if not voc_code in voc_list:
        new_voc = o.new_vocabulary(
            code=voc_code,
            description='TEST_VOC_DESCRIPTION',
            terms=[
                {"code": 'term_code1', "label": "term_label1",
                    "description": "term_description1"},
                {"code": 'term_code2', "label": "term_label2",
                    "description": "term_description2"},
                {"code": 'term_code3', "label": "term_label3",
                    "description": "term_description3"}
            ]
        )
        new_voc.save()
        print(f'Creating new vocabulary {new_voc.code}')
        created_voc = True
    else:
        new_voc = o.get_vocabulary(voc_code)
        print(f'Fetching exiting vocabulary {new_voc.code}')

    # CREATING NEW PROPERTY TYPES FOR TESTING WITH EVERY POSSIBLE DATA TYPE

    props_dict = {
        'TEST_BOOLEAN': 'BOOLEAN',
        'TEST_CONTROLLEDVOCABULARY': 'CONTROLLEDVOCABULARY',
        'TEST_HYPERLINK': 'HYPERLINK',
        'TEST_INTEGER': 'INTEGER',
        'TEST_MATERIAL': 'MATERIAL',
        'TEST_MULTILINE_VARCHAR': 'MULTILINE_VARCHAR',
        # 'TEST_OBJECT': 'OBJECT', CANNOT CREATE OBJECT PROPERTY TYPES FROM PYBIS Allowed values for enum dataType are: ['INTEGER', 'VARCHAR', 'MULTILINE_VARCHAR', 'REAL', 'TIMESTAMP', 'BOOLEAN', 'CONTROLLEDVOCABULARY', 'MATERIAL', 'HYPERLINK', 'XML']
        'TEST_REAL': 'REAL',
        'TEST_TIMESTAMP': 'TIMESTAMP',
        'TEST_VARCHAR': 'VARCHAR',
        'TEST_XML': 'XML',
    }
    pt_dict = {}
    pt_types = list(o.get_property_types().df['code'])
    created_properties = False

    for prop, val in props_dict.items():

        # print(f'Prop: {prop}')
        if not prop in pt_types:
            new_pt = o.new_property_type(
                code=prop,
                label=f'{prop.lower()}_label',
                description=f'{prop.lower()}_description',
                dataType=val,
                vocabulary=new_voc.code if prop == 'TEST_CONTROLLEDVOCABULARY' else None,
            )
            new_pt.save()
            print(f'Creating new property {new_pt.code}')
            created_properties = True
        else:
            new_pt = o.get_property_type(prop)
            print(f'Fetching existing property {new_pt.code}')

        pt_dict[new_pt.code] = new_pt

    # ASSIGNING THE NEWLY CREATED PROPERTIES TO THE NEW SAMPLE TYPE

    for i, p in enumerate(pt_dict.keys()):

        blockPrint()
        test_sample_type.assign_property(
            prop=p,
            section=f'Section_{p}',
            ordinal=(i+1),
            mandatory=True,
            initialValueForExistingEntities=f'Initial_Val_{p}',
            showInEditView=True,
            showRawValueInForms=True,
        )
        enablePrint()

        print(f'Assigned property {p} to {test_sample_type.code}')

    # TESTING DELETING SAMPLES

    if not created_sample:
        print(f'Deleting sample {test_sample_type}')
        test_sample_type.delete('testing deletion')

    # TESING DELETING PROPERTIES
    if not created_properties:

        for p_name, p in pt_dict.items():
            print(f'Deleting property {p_name}')
            p.delete('testing deletion')

    # TESTING DELETING VOCABULARIES
    if not created_voc:
        print(f'Deleting vocabulary {new_voc.code}')
        new_voc.delete('testing deletion')

    print('Done')


def main():
    o = ExpStep.connect_to_datastore()

    samples = o.get_samples(
        where={
            '$name': '3Dm3_0_1rpm_Vogel_2_7_T17_02'
        }
    )
    identifier = samples[0].identifier
    sample_dict = ExpStep.get_sample_dict(o, identifier)
    sample_list = list(o.get_sample_type(
        sample_dict['type']).get_property_assignments().df['propertyType'])

    testsample = ExpStep.load_sample(o, identifier)
    print(len(testsample.datasets))
    # print(testsample.get_property_types(o))


if __name__ == '__main__':
    main()
