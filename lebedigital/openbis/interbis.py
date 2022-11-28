import logging
import os
import sys
from getpass import getpass

import pandas as pd
from pybis import Openbis


class Interbis(Openbis):
    def __init__(self, url, verify_certificates=True, token=None, use_cache=True,
                 allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=False):
        super().__init__(url, verify_certificates, token, use_cache,
                         allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks)

    def connect_to_datastore(self, username: str = None, password: str = None):
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
                sys.exit(1)
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
        Generates a pandas Dataframe containing an import template for a given object type exiting in data store.
        Can be saved as an excel sheet, too, for manual data entry and then imported to an ExpStep object
            by ExpStep.import_from_template

        Args:
            sample_type (str): The type of the sample (ex. EXPERIMENTAL_STEP_SOMETHING)
            write (bool, optional): Specifies if the function should return the template as a pandas DataFrame or if it should be written to an excel file (True). Defaults to False.
            sheet_name (str, optional): The name of the sheet in an excel file. Defaults to 'metadata'.
            path (str, optional): Path where the excel sheet should be saved. Only use with write set to True. Defaults to ''.

        Raises:
            ValueError: Raises an error when write is set to True and no path is given

        Returns:
            pandas.DataFrame: returns a DataFrame with the import template
        """

        # Raises an error when no path given and write is set to True
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

        # Here we define some small internal functions to shorten the overall functions as these
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

    def get_sample_properties(self, sample_type: str) -> pd.DataFrame:
        """Returns a DataFrame of the sample properties with their descriptions, labels and other metadata

        Args:
            sample_type (str): The sample type for which properties should be fetched

        Returns:
            pd.DataFrame: DataFrame of all properties with their attributes
        """
        # Getting a list of all the sample's properties
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

        df = df.transpose().reset_index(drop=True)
        df = df.fillna(value="")
        return df

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
        """Returns the full identifier of the sample

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

    def get_collection_identifier(self, collection_code: str) -> str:
        """Returns the full identifier of the collection

        Args:
            collection_code (str): Name of the collection

        Returns:
            str: Identifier of the collection
        """
        collections_df = self.get_collections(props='$name').df
        collection_code = collection_code.upper()

        if len(collections_df.index):
            return collections_df[collections_df['identifier'].str.contains(
                collection_code)].iloc[0]['identifier']
        else:
            raise ValueError(f'No collection with name {collection_code} found')

    def exists_in_datastore(self, name: str) -> bool:
        """Checks if a sample with the given identifier exists in the openBIS datastore.

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
        """Used for automatically creating a sample type within the doit tasks.
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
            # print(f'Sample type {sample_code} already exists in the Datastore')
            new_sample_type = self.get_sample_type(sample_code)

        else:
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

        pt_dict = {}
        pt_types = list(self.get_property_types().df['code'])

        for prop, val in sample_properties.items():

            if not prop.upper() in pt_types:
                new_pt = self.new_property_type(
                    code=prop,
                    dataType=val[0],
                    label=val[1],
                    description=val[2],
                )
                new_pt.save()
            else:
                new_pt = self.get_property_type(prop)

            pt_dict[new_pt.code] = new_pt

        # ASSIGNING THE NEWLY CREATED PROPERTIES TO THE NEW SAMPLE TYPE

        for i, p in enumerate(pt_dict.keys()):
            block_print()
            new_sample_type.assign_property(
                prop=p,
                section='Metadata',
                ordinal=(i + 1),
                mandatory=True if p == '$NAME' else False,
                # initialValueForExistingEntities=f'Initial_Val_{p}',
                showInEditView=True,
                showRawValueInForms=True,
            )
            enable_print()

        logging.debug(f'Sample Type {sample_code} created.')
        return self.get_sample_type(sample_code)


def main():
    # I am just checking if the stuff i move still works, will delete when im done with porting
    o = Interbis('https://test.datastore.bam.de/openbis/')
    o.connect_to_datastore()

    # sample_type_name = 'EXPERIMENTAL_STEP_WILL_DELETE_LATER',

    # sample_type_dict = {
    #     'typeStr': ['VARCHAR', 'typeStr_label', 'typeStr_desc'],
    #     'typeInt': ['INTEGER', 'typeInt_label', 'typeInt_desc'],
    #     'typeFloat': ['REAL', 'typeFloat_label', 'typeFloat_desc'],
    #     'typeBoolean': ['BOOLEAN', 'typeBool_label', 'typeBool_desc'],
    # }

    # before_flag = False
    # sample_types_before = o.get_sample_types().df
    # print(sample_types_before)
    # if sample_type_name in sample_types_before.code.values:
    #     before_flag = True

    # sample_type_code = o.create_sample_type(
    #     sample_code='EXPERIMENTAL_STEP_WILL_DELETE_LATER',
    #     sample_prefix='ABCDEFG',
    #     sample_properties=sample_type_dict,
    # )

    # sample_type_code.delete(reason='aa')

    collection = o.get_collection_identifier('Lebedigital_Collection')
    print(collection)


if __name__ == '__main__':
    main()
    # full_emodul()
