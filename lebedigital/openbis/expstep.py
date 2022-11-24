import logging
import os
import shutil
import sys
from math import isnan
from pathlib import Path
from sys import exit

import pandas as pd
import yaml
from termcolor import colored

from lebedigital.openbis.interbis import Interbis
from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import \
    extract_metadata_emodulus


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
        """Creates an ExpStep Object

        Args:
            name (str, optional): name of the experiment. Defaults to ''.
            type (str, optional): type of the experiment, has to be a defined sample type in the Datastore. Defaults to ''.
            metadata (dict, optional): metadata of your data, is defined by the sample type. Defaults to None.
            space (str, optional): space in which the experiment should be saved. Defaults to "".
            project (str, optional): project in which the experiment should be saved. Defaults to "".
            collection (str, optional): collection in which the expeiment should be saved. Defaults to "".
            parents (list, optional): parents of the experiment. Defaults to None.
            identifier (str, optional): Verbose identifier of the sample, indicates the space and project. Defaults to ''.
            permId (str, optional): PermID identifier of the sample. Defaults to ''.
            sample_object (Pybis, optional): The pybis sample object of the sample. Defaults to None.
            datasets (list, optional): List of dataset objects saved under the experimental step. Defaults to None.
            dataset_codes (list, optional): The PermIDs of the datasets. Defaults to None.
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
        """Prints a description of the sample into the terminal/jupter output
        """
        print('\n')
        title = colored(self.name, 'cyan')
        title = '=============== SAMPLE "' + title + '" ================'
        print(title)
        print('\n')
        for param, value in vars(self).items():

            colparam = colored(
                param, 'green') if value else colored(param, 'red')

            if isinstance(value, str):
                print(f'{colparam: >25}: {value}' + '\n')
            elif isinstance(value, dict):
                print(f'{colparam: >25}:')
                for key, val in value.items():
                    colkey = colored(key, 'green', 'on_grey')
                    print(' '*20 + f'{colkey}: {val}')
                print('\n')
            else:
                print(f'{colparam: >25}: {value}' + '\n')

        print('=' * int(float(len(title)) * 0.9))

    def import_from_template(self, path_to_file: str):
        """Imports the data from the template into the samples metadata. To used with the generated import template

        Args:
            path_to_file (str): Path to the filled out import template
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
        types_df = o.get_sample_properties(self.type)
        types_dict = dict(zip(types_df.code, types_df.dataType))
        logging.debug(types_dict)
        types_dict = {k.lower(): conv_dict[v] for k, v in types_dict.items()}

        # Check if all keys are a subset of all possible keys
        logging.debug('Checking if keys are subset of defined parameters')
        if set(self.metadata.keys()).issubset(set(types_dict.keys())):
            for key, val in self.metadata.items():
                # Check if all values have the correct data type
                if not isinstance(val, types_dict[key]):
                    raise ValueError(
                        f'Type Checker: Type Mismatch. Entry "{key}" is of Type {type(val)} and should be {types_dict[key]}')
                # Check if float values are not NaN
                if isinstance(val, float) and isnan(val):
                    raise ValueError(
                        f'Type Checker: NaN values are not accepted by Openbis. Entry "{key} has a value "{val}".')
        else:
            # Keys in metadata and not in object definition -> raises an error for every such key
            for key in self.metadata.keys():
                if key.lower() not in types_dict.keys():
                    raise ValueError(
                        f'Type Checker: Unexpected parameter. Key "{key}" is not in {self.type} defined metadata params')

        # Checking if all values are defined
        needed_attrs = ['name', 'type', 'collection', 'space', 'project']
        for attr, value in self.__dict__.items():
            if attr in needed_attrs:
                if not value:
                    raise ValueError(
                        f'Type Checker: Undeclared Attribute. Attribute "{attr}" has not been delcared')

        # Printing warnings if parents are empty, does not break the function
        if not self.parents:
            logging.warning(
                'Type Checker: Parents attribute undeclared')

        # If you got to this line no gamebreaking bugs were found
        logging.info('Type Checker: Types are correctly set')
        
    @staticmethod
    def load_sample(o: Interbis, sample_identifier: str):
        """Loads an expstep from an experimental step in the datastore with its properties

        Args:
            o (Interbis): currently running openBIS interface
            sample_identifier (str): identifier of the sample

        Returns:
            ExpStep: ExpStep object containing metadata of the sample
        """

        # Getting the properties of the sample
        sample_dict = o.get_sample_dict(sample_identifier)

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
            identifier=sample_dict['identifier'],
            permId=sample_dict['permId'],
            sample_object=sample,
            datasets=sample_datasets,
            dataset_codes=sample_dataset_codes
        )

        return sample_step

    def fill_sample(self, o: Interbis, sample_identifier: str):
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

        return self

    def upload_expstep(self, o: Interbis, overwrite: bool = False) -> str:
        """Uploads the ExpStep object into the datastore

        Args:
            o (Interbis): Currently running openBIS interface
            name (str): Name of the sample
            overwrite (bool): Specifies if an existing Sample should be overwritten. Defaults to False.

        Raises:
            ValueError: ValueError if the Metadata is empty,
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

        # If a sample with the same name exists in the Datastore you fetch it instead of creating a new one
        if o.exists_in_datastore(self.name):
            print(f'Sample {self.name} already exists in Datastore')
            samples_df = o.get_samples(
                type=self.type,
                space=self.space,
                project=self.project,
                props="$name"
            ).df

            # Getting the identifier from the dataframe
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
                    type=self.type,
                    space=self.space,
                    collection=self.collection,
                    parents=self.parents,
                    # children=self.children,
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
                # code=self.code,
                type=self.type,
                space=self.space,
                collection=self.collection,
                parents=self.parents,
                # children=self.children,
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
        # Check if the sample exsts and delete it
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
        test = ''
        if '$name' in props:
            test = o.get_datasets(where={'$name': props['$name']})
        elif '$NAME' in props:
            test = o.get_datasets(where={'$NAME': props['$NAME']})
        else:
            raise KeyError('$name not defined in props')

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

    def download_datasets(self, o: Interbis, path: str, data_type: str = ''):
        """Downloads all datasets which are asigned to that sample

        Args:
            o (Interbis): Currently running openBIS interface
            path (str): Path where the datasets should be saved
            data_type (str): If specified will only download data sets of that type

        Raises:
            ValueError: Raises an error when no datasets are found under the sample
        """

        # If the sample has no datasets an error will be thrown
        if not len(self.datasets):
            raise ValueError('No Datasets found under the sample')

        file_plural = 'FILES' if len(self.datasets) > 1 else 'FILE-'

        print(
            f'----------DOWNLOADING {len(self.datasets)} {file_plural}----------\n')
        for dataset in self.datasets:
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

    def save_sample_yaml(self, yaml_path):
        """Saves a log file of a sample in openBIS. Used in doit upload

        Args:
            yaml_path (str): Path where the yaml file should be saved
        """

        modified_dict = self.__dict__

        # The python objects are ignored in the printout, too long and too useless
        modified_dict['sample_object'] = '--A PYTHON OBJECT WAS HERE--'
        modified_dict['datasets'] = '--A PYTHON OBJECT WAS HERE--'

        with open(yaml_path, 'w') as file:
            documents = yaml.dump(modified_dict, file)

# Here the tests are starting, can be run from the main function


def new_object_test():

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
        # 'TEST_MATERIAL': 'MATERIAL', WEIRD STUFF IS HAPPENNING WITH IT, WEIRD FORMATS WEIRD BEHAVIOUR ALL AROUND BAD STUFF
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


def download_datasets_test():

    def delete_folder(path):
        folder = path
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

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

    test_folder = '/home/ckujath/code/testing'
    testsample.download_datasets(o, test_folder)
    wait = input('Input something when ready to continue')
    delete_folder(test_folder)

    testsample.download_datasets(o, test_folder, 'ELN_PREVIEW')
    wait = input('Input something when ready to continue')
    delete_folder(test_folder)

    # print(testsample.get_property_types(o))


def import_template_test():

    o = ExpStep.connect_to_datastore()
    test_folder = '/home/ckujath/code/testing'

    ExpStep.gen_metadata_import_template(
        o, 'EXPERIMENTAL_STEP_TEST', True, 'metadata', f'{test_folder}/test_sheet.xlsx')


def load_yaml_test(path):

    o = ExpStep.connect_to_datastore()

    test_sample = ExpStep(
        name='test_sample',
    )

    test_sample.read_metadata_emodul(path)
    test_sample.sync_name('name')

    print(test_sample.__dict__)


def full_emodul():

    o = ExpStep.connect_to_datastore()

    metadata_path = '/home/ckujath/code/testing/Wolf 8.2 Probe 1.yaml'
    processed_data_path = '/home/ckujath/code/testing/Wolf 8.2 Probe 1.csv'
    preview_path = '/home/ckujath/code/testing/test_graph.png'

    emodul_sample = ExpStep(
        name='Wolf 8.2 Probe 1',
        space='CKUJATH',
        project='LEBEDIGITAL',
    )
    emodul_sample.read_metadata_emodul(metadata_path)

    emodul_sample.collection = o.find_collection(
        o, 'LEBEDIGITAL_COLLECTION', id_type=1)

    emodul_sample.sync_name(get_from='name')

    emodul_sample_type = ExpStep.create_sample_type_emodul(
        o,
        sample_code='EXPERIMENTAL_STEP_EMODUL',
        sample_prefix='EMODUL',
        sample_properties=emodul_sample.metadata,
    )

    emodul_sample.type = emodul_sample_type.code

    emodul_sample.upload_expstep(o)

    emodul_sample.upload_dataset(
        o,
        props={
            '$name': f'{emodul_sample.name}_processed',
            'files': processed_data_path,
            'data_type': 'PROCESSED_DATA'
        }
    )

    emodul_sample.upload_dataset(
        o,
        props={
            '$name': f'{emodul_sample.name}_preview',
            'files': preview_path,
            'data_type': 'PROCESSED_DATA'
        }
    )

    print(emodul_sample.info())


def responses_for_tests():
    o = ExpStep.connect_to_datastore()

    space = 'CKUJATH'
    project = 'TEST_AMBETON'
    collection = '/CKUJATH/TEST_AMBETON/VISKO_DATA_COLLECTION'

    # samples_df = o.get_samples(
    #     space=space,
    #     project=project,
    #     collection=collection,
    #     props=['$name']
    # ).df

    # samples_dict = pd.DataFrame.to_dict(samples_df)

    names = ExpStep.get_sample_names(o, space, project, collection)

    print(names)


def create_object_for_testing(space='CKUJATH', project='LEBEDIGITAL', collection='/CKUJATH/LEBEDIGITAL/LEBEDIGITAL_COLLECTION'):
    """Creates an object which can be used to test the functionality of the module

    Args:
        space (str, optional): Space NAME. Defaults to 'CKUJATH'.
        project (str, optional): Project NAME. Defaults to 'LEBEDIGITAL'.
        collection (str, optional): Collection IDENTIFIER. Defaults to '/CKUJATH/LEBEDIGITAL/LEBEDIGITAL_COLLECTION'.
    """

    o = ExpStep.connect_to_datastore()

    sample = ExpStep(
        name='testing_sample_lebedigital',
        type='EXPERIMENTAL_STEP_TEST',
        space=space,
        project=project,
        collection=collection,
    )

    sample.import_from_template('/home/ckujath/code/testing/test_sheet.xlsx')

    # need to update the type checker to auto convert from str "True" to bool True
    sample.metadata['test_boolean'] = True

    sample.upload_expstep(o)

    sample.info()

# Here the short tests end

def testing_sample_log():
    # setting up the test example
    input = '/home/ckujath/code/LebeDigital/tests/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.csv'
    mix_file = 'mix.dat'
    specimen_file = 'specimen.dat'

    target_data = {'experimentName': 'BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4',
                   'software_specification': 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A',
                   'operator_timestamp': '13:25:39',
                   'operator_date': '01.09.2014',
                   'tester_name': 'Kh',
                   'specimen_name': 'BA-Losert E-Modul 28d v. 04.08.14 Probe 4',
                   'remark': 'Kraftgeregelt 3,9 kN/s',
                   'weight': 5342.0,
                   'weight_unit': 'g',
                   'diameter': 98.6,
                   'length': 300.3,
                   'length_unit': 'mm',
                   'mix_file': '2014_08_05 Rezeptur_MI.xlsx'}

    # run extraction and getting a dictionary with metadata
    test_data = extract_metadata_emodulus(input, specimen_file, mix_file)
    yaml.dump(test_data, '/home/ckujath/code/testing/test_sheet.xlsx')


if __name__ == '__main__':
    # new_object_test()
    # download_datasets_test()
    # import_template_test()
    # load_yaml_test('/home/ckujath/code/testing/Wolf 8.2 Probe 1.yaml')
    # responses_for_tests()
    # full_emodul()
    # create_object_for_testing()
    testing_sample_log()
    print('Done')
