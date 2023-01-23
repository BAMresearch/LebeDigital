import yaml
import os
from pathlib import Path

conv_dict = {
    str: 'VARCHAR',
    float: 'REAL',
    int: 'INTEGER',
}

# parent directory of the minimum working example
this_file_path = Path(__file__)
ParentDir = os.path.dirname(this_file_path.parent.absolute())
mixture_output_directory = Path(ParentDir, 'mixture')
metadata_mixture_directory = Path(mixture_output_directory,
                                  'metadata_yaml_files')


def create_union_yaml(yaml_directory_path: Path, output_path: Path, mixture_code: str, defaults_dict: dict) -> dict:

    union_dict = defaults_dict

    for file in os.scandir(yaml_directory_path):
        with open(file, "r") as stream:
            current_dict = {key: val for key, val in dict(yaml.safe_load(stream)).items()
                            if key not in union_dict.keys()}
            current_dict = {f"{mixture_code}.{key}": [conv_dict[type(val)], key, key]
                            for key, val in current_dict.items()}
        union_dict = dict(current_dict, **union_dict)

    output_file = Path(output_path, "mixfile_union.yaml")

    with open(output_file, "w") as output:
        dump = yaml.dump(union_dict)
        output.write(dump)

    return union_dict


if __name__ == "__main__":
    create_union_yaml(metadata_mixture_directory, os.path.dirname(this_file_path))
