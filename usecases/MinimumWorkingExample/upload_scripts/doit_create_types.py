import logging
import os
from pathlib import Path
from typing import Union

import yaml
from pybis.entity_type import SampleType

from lebedigital.openbis.interbis import Interbis

conv_dict = {
    str: 'VARCHAR',
    float: 'REAL',
    int: 'INTEGER',
}

ingredient_code = "EMODUL_INGREDIENT"
ingredient_prefix = "EMODUL_ING"
ingredient_props = {"$name": ["VARCHAR", "Name", "Name"],
                    f"{ingredient_code}.bulk_density": ["REAL", "Bulk Density", "Bulk Density"],
                    f"{ingredient_code}.source": ["VARCHAR", "source", "source"]}


def create_required_sample_types(mixture_directory_path: Union[Path, str],
                                 emodul_directory_path: Union[Path, str],
                                 config: dict,
                                 default_props: dict,
                                 logging_path: Union[Path, str],
                                 ingredient_keywords: list):
    if config['runson'] == 'nodb':
        _create_logfiles(mixture_sample_type='RAN ON ACTIONS',
                         emodul_sample_type='RAN ON ACTIONS',
                         logging_path=logging_path)
        return

    # connecting to datastore
    o = Interbis(config['datastore_url'])
    o.connect_to_datastore(username=config['user'], password=config['pw'])

    ingredient_sample_type = o.create_sample_type(ingredient_code, ingredient_prefix, ingredient_props)

    # creating union dict
    mix_union_dict = _create_union_dict(
        yaml_directory_path=Path(mixture_directory_path),
        output_path=None,
        mixture_code=f"EXPERIMENTAL_STEP_{config['mixture_prefix']}",
        defaults_dict=default_props)

    # filtering out ingredient_keywords
    mix_union_dict = {key: val for key, val in mix_union_dict.items() if not [keyword for keyword in ingredient_keywords if keyword in key]}

    mixture_sample_type = o.create_sample_type(
        sample_code=f"EXPERIMENTAL_STEP_{config['mixture_prefix']}",
        sample_prefix=config['mixture_prefix'],
        sample_properties=mix_union_dict,
    )
    logging.debug(f'mixture type created: {mixture_sample_type.code}')

    metadata_path = Path(emodul_directory_path, os.fsdecode(os.listdir(emodul_directory_path)[0]))

    # CREATING EMODUL SAMPLE TYPE from single yaml (all equal)
    emodul_metadata_type_dict = _read_metadata_for_types(metadata_path, f"EXPERIMENTAL_STEP_{config['emodul_prefix']}",
                                                         default_props)

    # CREATING VOCABULARIES
    _create_unit_vocabularies(o)

    # Creating the emodul sample type with the formatted dict
    emodul_sample_type = o.create_sample_type(
        sample_code=f"EXPERIMENTAL_STEP_{config['emodul_prefix']}",
        sample_prefix=config['emodul_prefix'],
        sample_properties=emodul_metadata_type_dict,
    )
    logging.debug(f'emodul type created: {emodul_sample_type.code}')

    _create_logfiles(mixture_sample_type=mixture_sample_type.code,
                     emodul_sample_type=emodul_sample_type.code,
                     ingredient_sample_type=ingredient_sample_type.code,
                     logging_path=logging_path)

    o.logout()


def _create_union_dict(yaml_directory_path: Path, output_path: Union[Path, None], mixture_code: str,
                       defaults_dict: dict) -> dict:
    union_dict = defaults_dict
    for file in os.scandir(yaml_directory_path):
        with open(file, "r") as stream:
            current_dict = {key: val for key, val in dict(yaml.safe_load(stream)).items()
                            if key not in union_dict.keys()}
            current_dict = {f"{mixture_code}.{key}": [conv_dict[type(val)], key, key]
                            for key, val in current_dict.items()}
        union_dict = dict(current_dict, **union_dict)

    if output_path:
        output_file = Path(output_path, "mixfile_union.yaml")
        with open(output_file, "w") as output:
            output.write(yaml.dump(union_dict))

    return union_dict


def _read_metadata_for_types(yaml_path: Union[str, Path], sample_type_code: str, default_props: dict):
    """
    Reads the metadata from a single yaml file and converts the dict to a format accepted by the create sample type
    method

    Args:
        yaml_path: Path to single yaml file
        sample_type_code: Prefix for the unique (not default) metadata specified in this sample type
        default_props: Dictionary with default props which are shared between sample types

    Returns:

    """
    default_keys = default_props.keys()
    with open(yaml_path, 'r') as file:
        loaded = dict(yaml.safe_load(file))
        output_dict = default_props
        for key, val in loaded.items():
            if key not in default_keys:
                val = "" if val is None else val
                output_dict[f"{sample_type_code}.{key}".lower()] = [conv_dict[type(val)], key.split('.')[-1],
                                                                    key.split('.')[-1]]

        # Tweaking the dict to use our defined controlledvocabulary instead of varchar for units
        output_dict[f"{sample_type_code}.length_unit".lower()] = ['CONTROLLEDVOCABULARY', 'length_unit', 'length_unit',
                                                                  'LENGTH_VOCAB']
        output_dict[f"{sample_type_code}.weight_unit".lower()] = ['CONTROLLEDVOCABULARY', 'weight_unit', 'weight_unit',
                                                                  'WEIGHT_VOCAB']

        return output_dict


def _create_unit_vocabularies(o: Interbis, length_vocab_name: str = 'LENGTH_VOCAB',
                              weight_vocab_name: str = 'WEIGHT_VOCAB'):
    length_terms = [
        {"code": 'mm', "label": "mm", "description": "unit mm"},
        {"code": 'cm', "label": "cm", "description": "unit cm"},
        {"code": 'm', "label": "m", "description": "unit m"},
    ]
    weight_terms = [
        {"code": 'mg', "label": "mg", "description": "unit mg"},
        {"code": 'g', "label": "g", "description": "unit g"},
        {"code": 'kg', "label": "kg", "description": "unit kg"}
    ]

    try:
        o.get_vocabulary(length_vocab_name)
    except ValueError:
        # Possible issues when some different vocabulary already defined in database, the possible vocab values will not
        # be the same obv
        _ = o.new_vocabulary(
            code=length_vocab_name,
            terms=length_terms,
            description='Vocabulary for length units'
        ).save()

    try:
        o.get_vocabulary(weight_vocab_name)
    except ValueError:
        # Possible issues when some different vocabulary already defined in database, the possible vocab values will not
        # be the same obv
        _ = o.new_vocabulary(
            code=weight_vocab_name,
            terms=weight_terms,
            description='Vocabulary for weight units'
        ).save()


def _create_logfiles(mixture_sample_type: Union[SampleType, str], emodul_sample_type: Union[SampleType, str], ingredient_sample_type: Union[SampleType, str],
                     logging_path: Union[Path, str]):
    with open(Path(logging_path, 'mixture_sample_type.yaml'), 'w') as file:
        # file.write(yaml.dump({'Sample Type': mixture_sample_type}))
        print(mixture_sample_type, file=file)

    with open(Path(logging_path, 'emodul_sample_type.yaml'), 'w') as file:
        # file.write(yaml.dump({'Sample Type': emodul_sample_type}))
        print(emodul_sample_type, file=file)

    with open(Path(logging_path, 'ingredient_sample_type.yaml'), 'w') as file:
        # file.write(yaml.dump({'Sample Type': emodul_sample_type}))
        print(ingredient_sample_type, file=file)
