from lebedigital.openbis.interbis import Interbis
from pydantic import BaseModel, validator, create_model
from typing import Optional, Union
from pybis.entity_type import SampleType
from datetime import datetime
import pandas as pd

CONVERSION_DICT = {
    'BOOLEAN': bool,
    'CONTROLLEDVOCABULARY': None,
    'DATE': datetime,
    'TIMESTAMP': str,
    'VARCHAR': str,
    'INTEGER': int,
    'REAL': float
}


class AttrsValidatorr(BaseModel):
    space: str = ""
    project: str = ""
    experiment: str = ""
    parents: Optional[list[str]] = []
    children: Optional[list[str]] = []

    class Config:
        title = 'sampletypeschema'


def generate_validator(interbis_obj: Interbis, sample_type: Union[SampleType, str]):
    property_df = interbis_obj.get_sample_type_properties(sample_type)
    property_dict = pd.Series(property_df.dataType.values, index=property_df.code).to_dict()

    property_function_input = {key: (CONVERSION_DICT[val], ...) for key, val in property_dict.items()}

    return create_model(
        'Validator',
        **property_function_input
    )
