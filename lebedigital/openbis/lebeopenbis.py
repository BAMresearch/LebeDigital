from pybis import Openbis
from lebedigital.openbis.expstep import ExpStep
import pandas as pd
import os
from getpass import getpass
import logging
from pprint import pprint
from sys import exit

class LeBeOpenbis(Openbis):
    def __init__(self, url='https://test.datastore.bam.de/openbis/', verify_certificates=True, token=None, use_cache=True, allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=False):
        super().__init__(url, verify_certificates, token, use_cache, allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks)

    def connect_to_datastore(self, username: str = None, password: str = None):
        """
        Establishes a connection to an openBIS Datastore. If username/password are parsed then
        they take precedence over default logic.

        Args:
            username (str, optional): Username of an openBIS user. Defaults to None.
            password (str, optional): Password of an openBIS user. Defaults to None.
        """

        # If a session is already active (when running a notebook or script again) return the active object
        # Instead of connecting again

        if not self.is_session_active():
            logging.debug('Establishing a connection with ' + self.url)

            # We can parse the username and password as user input
            # If left empty the username will be grabbed from the os username
            # and password will have to be entered by the user
            if username:
                os.environ['OPENBIS_USERNAME'] = username
            else:
                os.environ['OPENBIS_USERNAME'] = os.getlogin()

            if password:
                os.environ['OPENBIS_PASSWORD'] = password
            else:
                os.environ['OPENBIS_PASSWORD'] = getpass("Give Password: ")

            try:
                self.login(os.environ['OPENBIS_USERNAME'], os.environ['OPENBIS_PASSWORD'])
            except ValueError:
                print("Wrong Credentials")
                exit(1)
        else:
            logging.debug('Connection already established with ' + self.url)

    def get_metadata_import_template(
        self,
        sample_type: str,
        write: bool = False,
        sheet_name: str = 'metadata',
        path: str = '',
    ) -> pd.DataFrame:
        """
        Generate a pandas Dataframe containing an import template for a given object type.
        Can be saved as an excel sheet to for manual data entry and then imported to an ExpStep object

        Args:
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

        meta_list = list(self.get_sample_type(
            sample_type).get_property_assignments().df['propertyType'])
        meta_dict = {idx: v for idx, v in enumerate(meta_list)}

        # Parse the dict into a dataframe

        meta_df = pd.DataFrame.from_dict(
            meta_dict, orient="index", columns=['Param'])

        # Get property descriptions
        props_list = list(self.get_sample_type(sample_type)
                          .get_property_assignments()
                          .df['propertyType'])

        # Build a dataframe with the properties and their descriptions
        df = pd.DataFrame()
        for prop in props_list:
            pt = self.get_property_type(prop)
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
            # If sheet exists remove it before writing the template
            if os.path.exists(path):
                os.remove(path)

            # Write the sheet
            with pd.ExcelWriter(path=path, engine='xlsxwriter') as writer:
                meta_df.to_excel(excel_writer=writer, sheet_name=sheet_name)

        else:
            return meta_df

    def get_overview(self, level: str, **kwargs) -> dict:
        """ Generates an overview for the samples stored in the datastore
            You need to provide the openBIS 'level' where you want the overview to start.
        Args:
            level (str): What entity should be the highest level of the overview (space/project/collection overview)

        Raises:
            ValueError: Raises an error when no correct value for the level was specified

        Returns:
            dict: Returns a dictionary with the overview
        """

        # We check if the openbis aliases were specified in the function arguments
        # if "experiment" in kwargs:
        #     kwargs["collection"] = kwargs["experiment"]
        #     kwargs.pop("experiment", None)

        # if "object" in kwargs:
        #     kwargs["sample"] = kwargs["object"]
        #     kwargs.pop("object", None)

        # Here we define internal functions to shorten the overall functions as these
        # are reused inside pretty often

        def get_space_names():
            space_names = self.get_spaces().df
            return list(space_names['code'].values)

        def get_project_names(space):
            project_names = self.get_projects(space=space).df
            return [name.split('/')[-1] for name in list(project_names['identifier'].values)]

        def get_collection_names(space, project):
            collection_names = self.get_experiments(space=space, project=project).df
            return [name.split('/')[-1] for name in list(collection_names['identifier'].values)]

        def get_sample_names_and_codes(space, project, collection):
            sample_names_df = self.get_samples(space=space, project=project, collection=collection, props=['$name']).df
            sample_names_df.columns = sample_names_df.columns.str.upper()
            sample_codes = [name.split('/')[-1] for name in list(sample_names_df['IDENTIFIER'].values)]
            sample_names = list(sample_names_df['$NAME'].values)
            return [f'{code} ({name})' for code, name in zip(sample_codes, sample_names)]

        # We go through all entries in for loops, the only difference between levels is where we start the loop
        # Watch out, takes ages to run because someone tested creating 500_000 pybis objects, ran 5 minutes for me
        if level == 'full':
            space_dict = {'DATASTORE': {}}

            # grabbing the names of all spaces
            space_names = get_space_names()

            for space in space_names:
                project_dict = {}

                # grabbing the names of all projects under space
                project_names = get_project_names(space)

                for project in project_names:
                    collection_dict = {}

                    # grabbing the names of all collections under project
                    collection_names = get_collection_names(space, project)

                    for collection in collection_names:
                        sample_list = []

                        # grabbing the names of all samples (with their codes) under collection
                        sample_names_and_codes = get_sample_names_and_codes(space, project, collection)
                        for sample in sample_names_and_codes:
                            sample_list.append(sample)

                        collection_dict[collection] = sample_list

                    project_dict[project] = collection_dict

                space_dict['DATASTORE'][space] = project_dict

            return space_dict

        elif level == 'space':
            space = kwargs.pop('space')
            project_dict = {space: {}}

            # grabbing the names of all projects under space
            project_names = get_project_names(space)

            for project in project_names:
                collection_dict = {}

                # grabbing the names of all collections under project
                collection_names = get_collection_names(space, project)

                for collection in collection_names:
                    sample_list = []

                    # grabbing the names of all samples (with their codes) under collection
                    sample_names_and_codes = get_sample_names_and_codes(space, project, collection)

                    for sample in sample_names_and_codes:
                        sample_list.append(sample)

                    collection_dict[collection] = sample_list

                project_dict[space][project] = collection_dict

            return project_dict

        elif level == 'project':
            space = kwargs.pop('space')
            project = kwargs.pop('project')
            collection_dict = {project: {}}

            # grabbing the names of all collections under project
            collection_names = get_collection_names(space, project)

            for collection in collection_names:
                sample_list = []

                # grabbing the names of all samples (with their codes) under collection
                sample_names_and_codes = get_sample_names_and_codes(space, project, collection)

                for sample in sample_names_and_codes:
                    sample_list.append(sample)

                collection_dict[project][collection] = sample_list

            return collection_dict

        elif level == 'collection':
            space = kwargs.pop('space')
            project = kwargs.pop('project')
            collection = kwargs.pop('collection')
            sample_list = []

            # grabbing the names of all samples (with their codes) under collection
            sample_names_and_codes = get_sample_names_and_codes(space, project, collection)

            for sample in sample_names_and_codes:
                sample_list.append(sample)

            sample_dict = {collection: sample_list}
            return sample_dict

        else:
            raise ValueError('No correct level specified')

    def get_sample_properties(self, sample_type: str) -> pd.DataFrame:
        """Returns a DataFrame of the sample properties with their descriptions, labels and other metadata

        Args:
            sample_type (str): The sample type for which properties should be fetched

        Returns:
            pd.DataFrame: DataFrame of all properties with their attributes
        """
        # Getting a list of all the samples properties
        props_list = list(self.get_sample_type(sample_type)
                          .get_property_assignments()
                          .df['propertyType'])

        df = pd.DataFrame()
        # Getting the metadata of every entry in props_list
        for prop in props_list:
            pt = self.get_property_type(prop)
            prop_series = pd.Series(pt.attrs.all())
            prop_df = pd.DataFrame(
                {prop: prop_series.values}, index=prop_series.index)

            # Combining the props together into a dataframe
            if df.empty:
                df = prop_df
            else:
                df = df.join(prop_df)

        return df.transpose()

    def get_sample_dict(self, identifier: str) -> dict:
        """Fetches a dictionary filled with information about the sample

        Very useful method actually if you already have the identifier. You can get all information about
        the sample from it like all metadata, collection, type, etc.

        Args:
            identifier (str): Identifier of the sample (identifier or permID)

        Returns:
            dict: dict containing all information about the sample
        """
        sample = self.get_sample(identifier)
        sample_info_dict = sample.attrs.all()
        sample_prop_dict = {key.upper(): val for key,
                            val in sample.p.all().items()}
        full_dict = (sample_info_dict | sample_prop_dict)

        # copying the information from key experiment to collection for ease of use
        full_dict['collection'] = full_dict.get('experiment')
        return full_dict
    
    def get_sample_identifier(self, name: str) -> str:
        sample_df = self.get_samples(where={'$name': name}).df
        if len(sample_df.index) > 1 or len(sample_df.index) == 0:
            raise ValueError(f'Couldnt find unique sample, the amount of samples with that name is {len(sample_df.index)}')
        
        return sample_df['identifier'].values[0]

    def exists_in_datastore(self, name: str) -> bool:
        """Checks wheter a sample with the given identifier exists in the openBIS datastore.

        Args:
            name (str): '$name' attribute of the sample

        Returns:
            bool: True when exists otherwise False
        """
        samples = self.get_samples(where={'$name': name}, props='$name').df
        
        # Warning if there are more than one entry in the Dataframe. Suggests something went wrong when uploading.
        df_length = len(samples.index)

        if not df_length:
            return False
        elif df_length > 1:
            logging.warning('More than one sample exists with the same name.')
            return True
        else:
            return True

    def upload_expstep(self, expstep: ExpStep, overwrite: bool = False) -> str:
        """
        Uploads the ExpStep object into the openBIS datastore

        Args:
            expstep (ExpStep): An ExpStep object which should be uploaded to the datastore
            overwrite (bool, optional): Specifies whether an exiting object should be overwritten. Defaults to False.

        Returns:
            str: Identifier of the uploaded sample
        """
        try:
            expstep.check_type(self)
        except Exception as e:
            logging.error(str(e))
            exit(1)

        # If a sample with the same name exists in the Datastore you fetch it instead of creating a new one
        if self.exists_in_datastore(expstep.name):
            logging.debug(f'Sample {expstep.name} already exists in Datastore')
            samples_df = self.get_samples(
                where={'$name': expstep.name},
                props="$name",
            ).df

            # Gettingthe identifier from the dataframe
            sample_identifier = samples_df['identifier'].values[0]

            # Overwriting the sample by deleting and uploading a new one
            # TODO: Overwriting samples by changing their metadata instead of deleting and uploading new sample

            if overwrite:
                logging.debug('Overwriting the sample')
                self.delete_expstep('overwriting sample')
                sample = self.new_sample(
                    type=expstep.type,
                    space=expstep.space,
                    collection=expstep.collection,
                    parents=expstep.parents,
                    props=expstep.metadata,
                )
                sample.save()
                sample_identifier = sample.identifier

        # No sample with the same name found -> uploading the sample
        else:
            logging.debug(f'Creating new sample {expstep.name}')
            sample = self.new_sample(
                type=expstep.type,
                space=expstep.space,
                collection=expstep.collection,
                parents=expstep.parents,
                props=expstep.metadata,
            )
            sample.save()
            sample_identifier = sample.identifier

        return sample_identifier

    def delete_expstep(self, identifier: str, reason: str):
        """Deletes a sample from the datastore

        Args:
            reason (str): Reason for deletion
        """
        # Check if the sample exsts and delete it
        try:
            self.get_sample(identifier).delete(reason)
        except Exception as e:
            logging.error(e)
            exit(1)

    def get_expstep(self, identifier: str) -> ExpStep:
        """Loads an expstep from an experimental step in the datastore with its properties

        Args:
            o (Openbis): currently running openbis instance
            sample_identifier (str): identifier of the sample

        Returns:
            ExpStep: ExpStep object containing metadata of the sample
        """

        # Getting the properties of the sample
        sample_dict = self.get_sample_dict(identifier)

        #  Getting a list of the properties only to flter them from the sample_dict
        props_list = list(self.get_sample_type(sample_dict['type'])
                          .get_property_assignments()
                          .df['propertyType'])

        # Getting the sample
        sample = self.get_sample(identifier, props='*')

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
        loaded_expstep = ExpStep(
            name=sample_dict['$NAME'],
            type=sample_dict['type'],
            space=sample_dict['identifier'].split('/')[1],
            project=sample_dict['identifier'].split('/')[2],
            collection=sample_collection,
            parents=sample_parents,
            metadata=sample_metadata,
            identifier=sample_dict['identifier'],
            permId=sample_dict['permId'],
            sample_object=sample,
            datasets=sample_datasets,
            dataset_codes=sample_dataset_codes
        )
        return loaded_expstep

    def upload_dataset(self, identifier: str, props: dict):
        """Uploads a dataset to the ExpStep.

        Requires a dictionary with $name, files and data_type

        Args:
            identifier (str): Identifier of the sample under which to store the dataset
            props (str, optional): Metadata of the dataset.

        Returns:
            str: Properties of the dataset
        """

        # Checking if the name of the dataset is included in props
        dataset_exists_test = ''
        if '$name' in props:
            dataset_exists_test = self.get_datasets(where={'$name': props['$name']})
        elif '$NAME' in props:
            dataset_exists_test = self.get_datasets(where={'$NAME': props['$NAME']})
        else:
            raise KeyError('$name not defined in props')

        # Checking if the files and data_type are specified in props
        if 'files' not in props:
            raise KeyError('files not specified in props')
        if 'data_type' not in props:
            raise KeyError('data_type not specified in props')

        files = props.pop('files')
        data_type = props.pop('data_type')

        # If a dataset with the same name was found in the datastore that dataset will be returned and none will be uploaded
        if dataset_exists_test:
            name = props['$name']
            print(
                f'Dataset(s) with the same name already present in the Datastore.\nTo upload the dataset you must first delete the other dataset with name {name}')
            ds = self.get_datasets(where={'$name': props['$name']})

        # Uploading the dataset
        else:
            ds = self.new_dataset(
                type=data_type,
                collection=self.get_sample_dict(identifier)['collection'],
                sample=identifier,
                files=files,
                props=props
            )
            ds.save()

        return ds

    def download_datasets(self, identifier: str, path: str, data_type: str = ''):
        """Downloads all datasets which are asigned to that sample

        Args:
            identifier (str): Identifier of the sample
            path (str): Path where the datasets should be saved
            data_type (str): If specified will only download data sets of that type

        Raises:
            ValueError: Raises an error when no datasets are found under the sample
        """

        # If the sample has no datasets an error will be thrown
        
        sample_object = self.get_sample(identifier)
        datasets = sample_object.get_datasets()
        
        if not len(datasets):
            raise ValueError('No Datasets found under the sample')

        file_plural = 'FILES' if len(datasets) > 1 else 'FILE-'

        print(
            f'----------DOWNLOADING {len(datasets)} {file_plural}----------\n')
        for dataset in datasets:
            # If data_type was specified download only the datastes with that data type
            if data_type:
                if dataset.type == data_type:
                    print(f'Downloading dataset {dataset.code}')
                    print(f'Files: {dataset.file_list}\n')
                    dataset.download(
                        destination=path,
                        create_default_folders=False,
                        wait_until_finished=False,
                    )

            # If data_type was NOT specfied download all datasets
            else:
                print(f'Downloading dataset {dataset.code}')
                print(f'Files: {dataset.file_list}\n')
                dataset.download(
                    destination=path,
                    create_default_folders=False,
                    wait_until_finished=False,
                )

        print('----------DOWNLOAD FINISHED----------')


def full_emodul():

    o = LeBeOpenbis()
    o.connect_to_datastore()

    metadata_path = '/home/ckujath/code/testing/Wolf 8.2 Probe 1.yaml'
    processed_data_path = '/home/ckujath/code/testing/Wolf 8.2 Probe 1.csv'
    preview_path = '/home/ckujath/code/testing/test_graph.png'

    emodul_sample = ExpStep(
        name='Wolf 8.2 Probe 1',
        space='CKUJATH',
        project='LEBEDIGITAL',
    )
    emodul_sample.collection = emodul_sample.find_collection(
        o, 'LEBEDIGITAL_COLLECTION', id_type=1)

    emodul_sample.sync_name(get_from='name')

    emodul_sample.read_metadata_emodul(metadata_path)

    emodul_sample_type = ExpStep.create_sample_type_emodul(
        o,
        sample_code='EXPERIMENTAL_STEP_EMODUL',
        sample_prefix='EMODUL',
        sample_properties=emodul_sample.metadata,
    )

    emodul_sample.type = emodul_sample_type.code

    emodul_sample.identifier = o.upload_expstep(emodul_sample)

    o.upload_dataset(
        emodul_sample.identifier,
        props={
            '$name': f'{emodul_sample.name}_processed',
            'files': processed_data_path,
            'data_type': 'PROCESSED_DATA'
        }
    )

    o.upload_dataset(
        emodul_sample.identifier,
        props={
            '$name': f'{emodul_sample.name}_preview',
            'files': preview_path,
            'data_type': 'PROCESSED_DATA'
        }
    )
    emodul_sample.datasets = o.get_datasets(emodul_sample.identifier)

    print(emodul_sample.info())

def main():
    # I am just checking if the stuff i move still works, will delete when im done with porting
    o = LeBeOpenbis()
    o.connect_to_datastore()

    identifier = '/CKUJATH/LEBEDIGITAL/EMODUL867319'
    sample_type = 'EXPERIMENTAL_STEP_EMODUL'

    df = o.get_metadata_import_template(sample_type)
    print(df)
    # overview = o.get_overview(level='project', space='CKUJATH', project='TEST_AMBETON')
    # pprint(overview)
    sample_props = o.get_sample_properties(sample_type)
    pprint(sample_props)
    pprint(o.get_sample_dict(identifier))

    expstep = o.get_expstep(identifier)
    expstep.info()
    
    visko_identifier = o.get_sample_identifier('3Dm3_0_1rpm_Vogel_2_7_T17_01')
    visko_step = o.get_expstep(visko_identifier)
    visko_step.info()
    
    o.download_datasets(identifier=visko_identifier, path='/home/ckujath/code/testing')

if __name__ == '__main__':
    # main()
    full_emodul()
