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


def create_union_yaml(yaml_directory_path: Path, output_path: Path) -> dict:
    union_dict = {"operator_date": ["DATE", "operator_date", "operator_date"],
                  "tester_name": ["VARCHAR", "tester_name", "tester_name"]}

    for file in os.scandir(yaml_directory_path):
        with open(file, "r") as stream:
            current_dict = yaml.safe_load(stream)
            current_dict = {f"EXPERIMENTAL_STEP_EMODUL_MIX.{key}": [conv_dict[type(val)], key, key]
                            if key not in union_dict.keys() else None
                            for key, val in current_dict.items()}
        union_dict = dict(current_dict, **union_dict)

    output_file = Path(output_path, "mixfile_union.yaml")

    with open(output_file, "w") as output:
        dump = yaml.dump(union_dict)
        output.write(dump)

    return union_dict


if __name__ == "__main__":
    create_union_yaml(metadata_mixture_directory, os.path.dirname(this_file_path))
