import logging
import os
import sys
from getpass import getpass
import json
import requests
import pandas as pd
from pybis import Openbis
from pybis.entity_type import SampleType
from typing import Optional, Union
from pydantic import create_model, AnyUrl, validator
from pydantic.main import ModelMetaclass
from enum import Enum
from dateutil.parser import parse


CONVERSION_DICT = {
    'BOOLEAN': bool,
    'DATE': str,
    'HYPERLINK': AnyUrl,
    'INTEGER': int,
    'MATERIAL': str,
    'MULTILINE_VARCHAR': str,
    'OBJECT': str,
    'REAL': float,
    'TIMESTAMP': str,
    'VARCHAR': str,
    'XML': str,  # TODO write a parser for XMLs
}


class Interbis(Openbis):
    def __init__(self, url, verify_certificates=False, token=None, use_cache=True,
                 allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=False):
        super().__init__(url, verify_certificates, token, use_cache,
                         allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks)

    def __dir__(self):
        return [
            "connect_to_datastore()",
            "get_metadata_import_template()",
            "get_overview()",
            "get_sample_type_properties()",
            "get_sample_dict()",
            "get_sample_identifier()",
            "get_dataset_identifier()",
            "get_collection_identifier()",
            "exists_in_datastore()",
            "create_sample_type()",
            "create_parent_hint()",
            "set_parent_annotation()",
            "get_parent_annotation()",
            "generate_typechecker()"
        ] + super().__dir__()

    def connect_to_datastore(self, username: Optional[str] = None, password: Optional[str] = None):
        """
        Establishes a connection to an openBIS Datastore. If username/password are parsed then
        they take precedence over default logic.

        Args:
            username (str, optional): Username of an openBIS user. Defaults to None.
            password (str, optional): Password of an openBIS user. Defaults to None.
        """

        # If a session is already active (when running a notebook or script again) return the active object
        # instead of connecting again

        if not self.is_session_active():
            logging.debug('Establishing a connection with ' + self.url)

            # We can parse the username and password as user input
            # If it is left empty the username will be grabbed from the os username
            # and the password will have to be entered by the user
            if username:
                os.environ['OPENBIS_USERNAME'] = username
            else:
                os.environ['OPENBIS_USERNAME'] = os.getlogin()

            if password:
                os.environ['OPENBIS_PASSWORD'] = password
            else:
                os.environ['OPENBIS_PASSWORD'] = getpass("Give Password: ")

            self.login(os.environ['OPENBIS_USERNAME'], os.environ['OPENBIS_PASSWORD'])

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
        Generates a pandas Dataframe containing an import template for a given object type existing in data store.
        Can be saved as an Excel sheet, too, for manual data entry and then imported to an ExpStep object
            by ExpStep.import_from_template

        Args:
            sample_type (str): The type of the sample (ex. EXPERIMENTAL_STEP_SOMETHING)
            write (bool, optional): Specifies if the function should return the template as a pandas DataFrame or if it
                should be written to an Excel file (True). Defaults to False.
            sheet_name (str, optional): The name of the sheet in an Excel file. Defaults to 'metadata'.
            path (str, optional): Path where the Excel sheet should be saved. Only use with write set to True.
                Defaults to ''.

        Raises:
            ValueError: Raises an error when write is set to True and no path is given

        Returns:
            pandas.DataFrame: returns a DataFrame with the import template
        """

        # Raises an error when no path is given and write is set to True
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

        # If write==True then write to file
        if write:
            # If sheet exists remove it before writing the template
            if os.path.exists(path):
                os.remove(path)

            # Write the sheet
            with pd.ExcelWriter(path=path, engine='xlsxwriter') as writer:
                meta_df.to_excel(excel_writer=writer, sheet_name=sheet_name)

        else:
            return meta_df

    def import_props_from_template(self, path_to_file: str) -> dict:

        df = pd.read_excel(path_to_file)
        df.columns = df.columns.str.lower()

        # if some params were left out then the row gets dropped
        df = df.dropna()

        # not useful for uploading to the datastore only for the user to know what to input
        df.pop('label')
        df.pop('description')

        # turn the df columns into a dict
        metadata = dict(zip([par.lower() for par in df.param], df.value))

        for key in metadata.keys():
            property_type = self.get_property_type(key)
            property_data_type = property_type.dataType
            if property_data_type == "BOOLEAN":
                if metadata[key] in ["True", "False"]:
                    metadata[key] = eval(metadata[key])
                else:
                    raise KeyError(f"Invalid key \"{key}\" for dataType BOOLEAN. "
                                   f"Expected \"True\" or \"False\", found {metadata[key]}")
            elif property_data_type == "INTEGER":
                metadata[key] = int(metadata[key])

        return metadata

    def get_overview(self, level: str, **kwargs) -> dict:
        """
        Generates an overview for the samples stored in the datastore
        You need to provide the openBIS 'level' where you want the overview to start.

        Args:
            level (str): What entity should be the highest level of the overview (space/project/collection overview)

        Raises:
            ValueError: Raises an error when no correct value for the level was specified

        Returns:
            dict: Returns a dictionary with the overview
        """

        # Here we define some small internal functions to shorten the overall functions as these
        # are reused inside pretty often

        def get_space_names():
            space_names = self.get_spaces().df
            return list(space_names['code'].values)

        def get_project_names(space):
            project_names = self.get_projects(space=space).df
            return list(project_names['code'].values)

        def get_collection_names(space, project):
            collection_names = self.get_experiments(space=space, project=project).df
            return [name.split('/')[-1] for name in list(collection_names['identifier'].values)]

        def get_sample_names_and_codes(space, project, collection):
            sample_names_df = self.get_samples(space=space, project=project, collection=collection, props=['$name']).df
            sample_names_df.columns = sample_names_df.columns.str.upper()
            sample_codes = [name.split('/')[-1] for name in list(sample_names_df['IDENTIFIER'].values)]
            sample_names = list(sample_names_df['$NAME'].values)
            return [f'{code} ({name})' for code, name in zip(sample_codes, sample_names)]

        # We go through all entries in loops, the only difference between levels is where we start the loop
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

    def get_sample_type_properties(self, sample_type: Union[SampleType, str]) -> pd.DataFrame:
        """
        Returns a DataFrame of the sample properties with their descriptions, labels and other metadata

        Args:
            sample_type (str): The sample type for which properties should be fetched

        Returns:
            pd.DataFrame: DataFrame of all properties with their attributes
        """
        if isinstance(sample_type, str):
            sample_type = self.get_sample_type(sample_type)

        # Getting a list of all the sample's properties
        props_list = list(sample_type
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

        df = df.transpose().reset_index(drop=True)
        df = df.fillna(value="")
        return df

    def get_sample_dict(self, identifier: str) -> dict:
        """
        Fetches a dictionary filled with information about the sample

        Very useful method actually if you already have the identifier. You can get all information about
        the sample from it like all metadata, collection, type, etc.

        Args:
            identifier (str): Identifier of the sample (identifier or permID)

        Returns:
            dict: dict containing all information about the sample
        """
        sample = self.get_sample(identifier)
        sample_info_dict = sample.attrs.all()
        sample_prop_dict = {key.upper(): val for key, val in sample.p.all().items()}
        full_dict = (sample_info_dict | sample_prop_dict)

        # copying the information from key experiment to collection for ease of use
        full_dict['collection'] = full_dict.get('experiment')
        return full_dict

    def get_sample_identifier(self, name: str) -> str:
        """
        Returns the full identifier of the sample

        Args:
            name (str): '$name' attribute of the sample

        Returns:
            str: Identifier of sample
        """
        sample_df = self.get_samples(where={'$name': name}).df
        if len(sample_df.index) > 1 or len(sample_df.index) == 0:
            raise ValueError(
                f'Could not find unique sample, the amount of samples with that name is {len(sample_df.index)}')

        return sample_df['identifier'].values[0]

    def get_dataset_permid(self, name: str) -> list[str]:
        """
        Returns the permId of the dataset

        Args:
            name (str): '$name' attribute of the dataset

        Returns:
            list[str]: List of the dataset identifiers
        """
        df = self.get_datasets(**{"$NAME": name}).df
        if len(df) == 0:
            logging.warning("No dataset found")
            return []
        elif len(df) > 1:
            logging.warning(f"More than one dataset with the name {name} found")

        return df['permId'].to_list()

    def get_collection_identifier(self, collection_code: str) -> str:
        """
        Returns the full identifier of the collection provided only one collection
        with the specified name exists

        Args:
            collection_code (str): Name of the collection

        Returns:
            str: Identifier of the collection
        """
        collections_df = self.get_collections().df
        collections_df['code'] = collections_df['identifier'].str.split("/").str[-1]

        collection_code = collection_code.upper()

        relevant_rows = collections_df.loc[collections_df['code'] == collection_code]

        if relevant_rows.empty:
            raise ValueError(f'No collection with name {collection_code} found')

        return relevant_rows['identifier'].values[0]

    def exists_in_datastore(self, name: str) -> bool:
        """
        Checks if a sample with the given identifier exists in the openBIS datastore.

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

    def create_sample_type(self, sample_code: str, sample_prefix: str, sample_properties: dict):
        """
        Used for automatically creating a sample type within the doit tasks.
        May be useful for automating upload, less customisation options than creating them "by hand" though

        Args:
            sample_code (str): The code of the sample is the name of the sample, for example EXPERIMENTAL_STEP_EMODUL
            sample_prefix (str): The prefix of the new sample,
                will appear before the code of the created sample as in CODE12345
            sample_properties (dict): The object properties which should be assigned to the sample.
                If the property code is not found in the datastore a new property will be created.

        Returns:
            Pybis: returns the pybis sample type object
        """

        # SUPPRESSING PRINTS FOR ASSIGNING PROPERTIES
        # Disable
        def block_print():
            sys.stdout = open(os.devnull, 'w')

        # Restore
        def enable_print():
            sys.stdout = sys.__stdout__

        sample_types = self.get_sample_types().df
        if sample_code in list(sample_types['code']):
            new_sample_type = self.get_sample_type(sample_code)

        else:
            # Create the SampleType object
            new_sample_type = self.new_sample_type(
                code=sample_code,
                generatedCodePrefix=sample_prefix,
                # description='Testing Experimental Step', DOES NOT WORK WITH PYBIS, PYBIS DOES NOT ACCEPT
                # DESCRIPTION ARGUMENT
                autoGeneratedCode=True,
                subcodeUnique=False,
                listable=True,
                showContainer=False,
                showParents=True,
                showParentMetadata=True,
                validationPlugin='EXPERIMENTAL_STEP.date_range_validation'
            )
            new_sample_type.save()

        pt_dict = self.create_property_types(sample_properties)

        # ASSIGNING THE NEWLY CREATED PROPERTIES TO THE NEW SAMPLE TYPE

        for i, p in enumerate(pt_dict.keys(), start=1):
            block_print()
            new_sample_type.assign_property(
                prop=p,
                section='Metadata',
                ordinal=i,
                mandatory=False,
                # initialValueForExistingEntities=f'Initial_Val_{p}',
                showInEditView=True,
                showRawValueInForms=True,
            )
            enable_print()

        logging.debug(f'Sample Type {sample_code} created.')
        return self.get_sample_type(sample_code)

    def create_property_types(self, sample_properties: dict) -> dict:

        # Create a dictionary with the samples
        pt_dict = {}

        # Get all possible property types from the datastore
        pt_types = list(self.get_property_types().df['code'])

        for prop, val in sample_properties.items():

            if not prop.upper() in pt_types:
                logging.debug(f"Creating new property type {prop.upper()}")
                new_pt = self.new_property_type(
                    code=prop,
                    dataType=val[0],
                    label=val[1],
                    description=val[2],
                    # The vocabulary needs to exist at this point
                    vocabulary=val[3] if val[0] == 'CONTROLLEDVOCABULARY' else None
                )
                new_pt.save()
            else:
                logging.debug(f"Fetching existing property type {prop.upper()}")
                new_pt = self.get_property_type(prop)

            pt_dict[new_pt.code] = new_pt

        return pt_dict

    def create_parent_hint(
        self,
        sample_type: Union[str, SampleType],
        label: str,
        parent_type: Union[str, SampleType],
        min_count: Optional[int] = None,
        max_count: Optional[int] = None,
        annotation_properties: Optional[list] = None
    ):
        """
        Method for creating parent hints with comments, has to be set before the parent annotation can be specified, similar to SampleType before Sample
        """

        if isinstance(sample_type, SampleType):
            sample_type = sample_type.code

        if isinstance(parent_type, SampleType):
            parent_type = parent_type.code

        settings_sample = self.get_sample("/ELN_SETTINGS/GENERAL_ELN_SETTINGS")
        settings = json.loads(settings_sample.props["$eln_settings"])

        hint = {
            "LABEL": label,
            "TYPE": parent_type,
        }

        if min_count:
            assert min_count >= 0, "min_count can not be negative"
            hint["MIN_COUNT"] = min_count

        if max_count:
            assert max_count >= 0, "max_count can not be negative"
            hint["MAX_COUNT"] = max_count

        if annotation_properties:
            hint["ANNOTATION_PROPERTIES"] = annotation_properties
        else:
            hint["ANNOTATION_PROPERTIES"] = [{
                "TYPE": "ANNOTATION.SYSTEM.COMMENTS",
                "MANDATORY": False,
            }]

        # For the case that no settings for the sample type have been set before
        settings["sampleTypeDefinitionsExtension"].setdefault(sample_type, {})
        settings["sampleTypeDefinitionsExtension"][sample_type].setdefault("SAMPLE_PARENTS_HINT", [])

        if hint not in settings["sampleTypeDefinitionsExtension"][sample_type]["SAMPLE_PARENTS_HINT"]:
            settings["sampleTypeDefinitionsExtension"][sample_type]["SAMPLE_PARENTS_HINT"].append(hint)

        settings_sample.props["$eln_settings"] = json.dumps(settings)

        settings_sample.save()

    def get_parent_annotations(self, sample_identifier: str) -> dict:
        sample = self.get_sample(sample_identifier)
        return sample.data["parentsRelationships"]

    def set_parent_annotation(self, child_sample: str, parent_sample: str, comment: str, set_property: Optional[str] = None):
        """
        Sets the ANNOTATION.SYSTEM.COMMENTS field for an existing parent-child relationship

        Args:
            child_sample(str): The identifier of the child sample in the relationship
            parent_sample(str): The identifier of the parent sample in the relationship
            comment(str): The value of the annotation
            set_property(str): The property of the annotation to be commented
        """
        def combine_urls(base_url, relative_url):
            # Remove any trailing or leading slashes
            base_url = base_url.strip('/')
            relative_url = relative_url.strip('/')

            # Remove the last object in the base URL
            base_url_parts = base_url.split('/')
            base_url_parts.pop()
            base_url = '/'.join(base_url_parts)

            # Combine the two URLs
            full_url = f'{base_url}/{relative_url}'

            return full_url

        child_sample_permid = self.get_sample(child_sample).permId

        if not set_property:
            set_property = "ANNOTATION.SYSTEM.COMMENTS"

        request = {'method': 'updateSamples',
                   'params': [self.token,
                              [{'@id': 0,
                                'properties': {},
                                  'relationships': {parent_sample: {'@id': 1,
                                                                    'parentAnnotations': {'@id': 2,
                                                                                          'actions': [{'@id': 3,
                                                                                                       'items': [{set_property: str(comment)}],
                                                                                                       '@type': 'as.dto.common.update.ListUpdateActionAdd'}],
                                                                                          '@type': 'as.dto.common.update.ListUpdateMapValues'},
                                                                    '@type': 'as.dto.common.update.RelationshipUpdate'}},
                                  'sampleId': {'@id': 4,
                                               'permId': child_sample_permid,
                                               '@type': 'as.dto.sample.id.SamplePermId'},
                                  '@type': 'as.dto.sample.update.SampleUpdate'}]],
                   'id': '1',
                   'jsonrpc': '2.0'}

        response = requests.post(url=combine_urls(self.url, self.as_v3), json=request, verify=self.verify_certificates).json()
        return response

    def _get_datatype_conversion(self, property_name: str, property_datatype: str):
        """
        Converts the openbis datatypes into python datatypes if possible, else a custom datatype
        Use with the `generate_typechecker` method
        """
        # if the prop is in the dict then it is not a CONTROLLED_VOCABULARY
        # the CONTROLLED_VOCABULARY property type is the only one which cant be simply mapped to an exisiting
        # python data type as we are handling it using an Enum which we create below
        if property_datatype in CONVERSION_DICT:
            return CONVERSION_DICT[property_datatype]

        vocabulary_name = self.get_property_type(property_name).vocabulary
        vocabulary_df = self.get_vocabulary(vocabulary_name).get_terms().df
        vocabulary_term_list = vocabulary_df['code'].to_list()
        vocabulary_enum_dict = dict(zip(vocabulary_term_list, vocabulary_term_list))

        return Enum('Vocabulary', vocabulary_enum_dict)

    def generate_typechecker(self, sample_type: Union[SampleType, str]) -> ModelMetaclass:
        """
        Generates a pydantic-based typechecker with property types saved in openbis for a given sample.

        Args:
            sample_type(str | SampleType): The code of the sample the typechecker should be created for
        Returns:
            ModelMetaclass: A dynamically created pydantic model
        """

        if not isinstance(sample_type, SampleType):
            sample_type = self.get_sample_type(sample_type)

        mandatory_props_df = sample_type.get_property_assignments().df
        mandatory_props_dict = pd.Series(mandatory_props_df.mandatory.values, index=mandatory_props_df.propertyType).to_dict()
        mandatory_props = [key.lower() for key, val in mandatory_props_dict.items() if val]

        property_df = self.get_sample_type_properties(sample_type)
        property_dict = pd.Series(property_df.dataType.values, index=property_df.code).to_dict()
        property_dict = {key.lower(): val for key, val in property_dict.items()}

        property_function_input = {key: (self._get_datatype_conversion(key, val), None if key not in mandatory_props else ...) for key, val in property_dict.items()}

        datetime_props = {key: val for key, val in property_dict.items() if val == "DATE" or val == "TIMESTAMP"}
        controlledvocabulary_props = {key: val for key, val in property_dict.items() if val == "CONTROLLEDVOCABULARY"}

        validators = {}

        if datetime_props:

            def date_correct_format(cls, v):
                return parse(v).strftime("%Y-%m-%d")

            def timestamp_correct_format(cls, v):
                return parse(v).strftime("%Y-%m-%d %H:%M")

            validators = validators | {f"{key}_validator": validator(key, pre=True, allow_reuse=True)(date_correct_format if val == 'DATE' else timestamp_correct_format) for key, val in datetime_props.items()}

        if controlledvocabulary_props:

            def props_to_uppercase(cls, v):
                return str(v).upper()

            validators = validators | {f"{key}_validator": validator(key, pre=True, allow_reuse=True)(props_to_uppercase) for key, val in controlledvocabulary_props.items()}

        class Config:
            extra = "forbid"
            use_enum_values = True

        return create_model(
            'SampleType_Props_Validator',
            **property_function_input,
            __config__=Config,
            __validators__=validators
        )
