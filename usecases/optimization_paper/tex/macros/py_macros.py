


stream = file('py_macros.yaml', 'r')
dictionary = yaml.load_all(stream)

for key in dict:
        print key, dict[key]
