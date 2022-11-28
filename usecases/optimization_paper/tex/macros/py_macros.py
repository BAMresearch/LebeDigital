import yaml

def py_macros(file_name):
    file_name = str(file_name) # convert pathlib object to useful string

    # read dictionaries
    with open(file_name + '.yaml') as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)

    # write tex macros
    with open(file_name + '.tex', 'w+') as f:
            for dictionary in data:
                for key in data[dictionary]:
                    f.write(f'\\newcommand{{\\{key}}}{{{data[dictionary][key]}}}\n')


if __name__ == "__main__":
    py_macros('py_macros')