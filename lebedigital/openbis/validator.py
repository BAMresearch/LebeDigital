from lebedigital.openbis.interbis import Interbis
from pydantic import BaseModel, create_model, Field
from typing import Optional, Union
from pybis.entity_type import SampleType
from datetime import datetime
import os
import pandas as pd
import pprint

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

    name_prop = property_dict.pop('$NAME')

    property_function_input = {key: (CONVERSION_DICT[val], None) for key, val in property_dict.items()}

    property_function_input['$name'] = (CONVERSION_DICT[name_prop], ...)

    pp = pprint.PrettyPrinter(depth=4)
    pp.pprint(property_function_input)

    class Config:
        extra = "forbid"

    return create_model(
        'Validator',
        **property_function_input,
        __config__=Config
    )


if __name__ == "__main__":
    o = Interbis(os.environ['OPENBIS_URL'])
    # print(o.is_session_active())

    new_sample = o.new_sample('EXPERIMENTAL_STEP_EMODUL_MIX')

    new_sample_props = {
        '$name': 'test_name',
        'experimental_step_emodul_mix.water_cement_ratio': 1.23,
        'experimental_step_emodul_mix.admixture--volume': 'ABC',
    }

    new_sample.props = new_sample_props

    Model = generate_validator(o, 'EXPERIMENTAL_STEP_EMODUL_MIX')
    print(Model)
    print(Model.schema_json(indent=2))

    Model.validate(new_sample.props.all_nonempty())
