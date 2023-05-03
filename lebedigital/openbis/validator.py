from lebedigital.openbis.interbis import Interbis
from pydantic import BaseModel, create_model, Field, AnyUrl
from typing import Optional, Union, List
from enum import Enum
from pybis.entity_type import SampleType
from datetime import datetime
import os
import pandas as pd
import pprint

CONVERSION_DICT = {
    'BOOLEAN': bool,
    'DATE': datetime,
    'HYPERLINK': AnyUrl,
    'INTEGER': int,
    'MATERIAL': str,
    'MULTILINE_VARCHAR': None,  # TODO find out how multilines are saved
    'OBJECT': str,
    'REAL': float,
    'TIMESTAMP': str,
    'VARCHAR': str,
    'XML': None  # TODO find out how xmls are saved
}


def get_conversion(interbis_obj: Interbis, property_name: str, property_datatype: str):
    # if the prop is in the dict then it is not a CONTROLLED_VOCABULARY
    if property_datatype in CONVERSION_DICT:
        return CONVERSION_DICT[property_datatype]

    vocabulary_name = o.get_property_type(property_name).vocabulary
    vocabulary_df = o.get_vocabulary(vocabulary_name).get_terms().df
    vocabulary_term_list = vocabulary_df['code'].to_list()
    vocabulary_enum_dict = dict(zip(vocabulary_term_list, vocabulary_term_list))

    return list[Enum('Vocabulary', vocabulary_enum_dict)]


def generate_validator(interbis_obj: Interbis, sample_type: Union[SampleType, str]):
    property_df = interbis_obj.get_sample_type_properties(sample_type)
    property_dict = pd.Series(property_df.dataType.values, index=property_df.code).to_dict()

    name_prop = property_dict.pop('$NAME')

    property_function_input = {key.lower(): (get_conversion(o, key, val), None) for key, val in property_dict.items()}

    property_function_input['$name'] = (get_conversion(o, '$name', name_prop), ...)

    pp = pprint.PrettyPrinter(depth=4)
    pp.pprint(property_function_input)

    class Config:
        extra = "forbid"

    return create_model(
        'SampleType_Props_Validator',
        **property_function_input,
        __config__=Config
    )


if __name__ == "__main__":
    o = Interbis(os.environ['OPENBIS_URL'])
    # print(o.is_session_active())

    new_sample = o.new_sample('EXPERIMENTAL_STEP_EMODUL')

    new_sample_props = {
        '$name': 'test_name',
        'experimental_step_emodul.weight_unit': 'GIGAGRAM',
        'experimental_step_emodul.length_unit': 'M',
    }

    new_sample.set_props(new_sample_props)

    Model = generate_validator(o, 'EXPERIMENTAL_STEP_EMODUL')
    print(Model)
    # print(Model.schema_json(indent=2))

    Model.validate(new_sample.props.all_nonempty())
