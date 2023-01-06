import logging
from math import isnan
from sys import exit
from typing import Iterable

import pandas as pd
import yaml
from termcolor import colored

from lebedigital.openbis.interbis import Interbis
from pybis.sample import Sample


class ISample(Sample):

    def import_from_template(self, path_to_file: str):
        """
        Imports the data from the template into the sample's metadata. To used with the generated import template
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
        metadata = dict(zip([par.lower() for par in df.param], df.value))

        for key, val in metadata:
            setattr(self.p, key, val)

    def check_type(self, o: Interbis):
        """
        Dry checking if the types of metadata entries are correct and if the upload step will run as expected.

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
