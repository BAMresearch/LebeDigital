import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Union

import yaml
from pybis.entity_type import SampleType
from upload_scripts.doit_upload import _read_metadata

from lebedigital.openbis.interbis import Interbis

conv_dict = {
    str: 'VARCHAR',
    float: 'REAL',
    int: 'INTEGER',
}


def create_required_sample_types(mixture_directory_path: Union[Path, str],
                                 emodul_directory_path: Union[Path, str],
                                 config: dict,
                                 default_props: dict,
                                 logging_path: Union[Path, str]):
    if config['runson'] == 'actions':
        _create_logfiles(mixture_sample_type='RAN ON ACTIONS',
                         emodul_sample_type='RAN ON ACTIONS',
                         logging_path=logging_path)
        return

    # connecting to datastore
    o = Interbis(config['datastore_url'])
    o.connect_to_datastore(username=config['user'], password=config['pw'])

    # creating union dict
    mix_union_dict = _create_union_yaml(
        yaml_directory_path=Path(mixture_directory_path),
        output_path=None,
        mixture_code=f"EXPERIMENTAL_STEP_{config['mixture_prefix']}",
        defaults_dict=default_props)

    mixture_sample_type = o.create_sample_type(
        sample_code=f"EXPERIMENTAL_STEP_{config['mixture_prefix']}",
        sample_prefix=config['mixture_prefix'],
        sample_properties=mix_union_dict,
    )
    logging.debug(f'mixture type created: {mixture_sample_type.code}')

    metadata_path = Path(emodul_directory_path, os.fsdecode(os.listdir(emodul_directory_path)[0]))

    # CREATING EMODUL SAMPLE TYPE from single yaml (all equal)
    emodul_metadata = _read_metadata(metadata_path, f"EXPERIMENTAL_STEP_{config['emodul_prefix']}", default_props)
    emodul_metadata_type_dict = _reformat_sample_dict(loaded_dict=emodul_metadata, defaults_dict=default_props)

    # Creating the emodul sample type with the formatted dict
    emodul_sample_type = o.create_sample_type(
        sample_code=f"EXPERIMENTAL_STEP_{config['emodul_prefix']}",
        sample_prefix=config['emodul_prefix'],
        sample_properties=emodul_metadata_type_dict,
    )
    logging.debug(f'emodul type created: {emodul_sample_type.code}')

    _create_logfiles(mixture_sample_type=mixture_sample_type.code,
                     emodul_sample_type=emodul_sample_type.code,
                     logging_path=logging_path)

    o.logout()


def _create_union_yaml(yaml_directory_path: Path, output_path: Union[Path, None], mixture_code: str,
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


def _reformat_sample_dict(loaded_dict: dict, defaults_dict: dict):
    conv_dict = {
        str: 'VARCHAR',
        float: 'REAL',
        int: 'INTEGER',
    }

    #output_dict = defaultdict(lambda: ["NA", "NA", "NA"])
    #output_dict['$name'] = ['VARCHAR', 'Name', 'Name']
    output_dict = defaults_dict # default properties

    for key, val in loaded_dict.items():
        val = "" if val is None else val
        if key not in output_dict.keys():
            output_dict[key.lower()] = [conv_dict[type(val)], key.split('.')[-1], key.split('.')[-1]]

    print('CHECK',output_dict)
    # input()
    #adding a type for mix_file
    output_dict['experimental_step_emodul.mix_file']=['VARCHAR','mix_file','mix_file']
    return dict(output_dict)


def _create_logfiles(mixture_sample_type: Union[SampleType, str], emodul_sample_type: Union[SampleType, str],
                     logging_path: Union[Path, str]):
    with open(Path(logging_path, 'mixture_sample_type.yaml'), 'w') as file:
        # file.write(yaml.dump({'Sample Type': mixture_sample_type}))
        print(mixture_sample_type, file=file)

    with open(Path(logging_path, 'emodul_sample_type.yaml'), 'w') as file:
        # file.write(yaml.dump({'Sample Type': emodul_sample_type}))
        print(emodul_sample_type, file=file)
