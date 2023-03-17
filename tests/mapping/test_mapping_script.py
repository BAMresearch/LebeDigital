import unittest
import rdflib
import datetime
import os
import yaml
from pathlib import Path

from lebedigital.mapping.emodul_mapping import generate_knowledge_graph, get_date_time_value


def query_and_test_result(typeOfQuery, predicate, prop, g):
    print(f"""Query for {typeOfQuery}:""")
    q = f"""
            SELECT ?o WHERE {{
                    {predicate} {prop} ?o .
            }}
            """
    print(q)
    result = g.query(q)
    return result


def print_next():
    print("Pass")
    print("###################Next Query#####################################")


def print_first():
    print("###################First Query####################################")


def test_generate_knowledge_graph():
    metadataPath = 'testMetaData.yaml'
    filename = "knowledgeGraph.ttl"
    print("Start testing:")
    # Generate a testMetadata file
    print("Generate metadata file for test:")
    target_data = {'experimentName': 'BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4',
                   'software_specification': 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A',
                   'operator_timestamp': '13:25:39',
                   'operator_date': '01.09.2014',
                   'tester_name': 'Kh',
                   'specimen_name': 'BA-Losert E-Modul 28d v. 04.08.14 Probe 4',
                   'remark': 'Kraftgeregelt 3,9 kN/s',
                   'weight': 5342.0,
                   'weight_unit': 'g',
                   'diameter': 98.6,
                   'length': 300.3,
                   'length_unit': 'mm',
                   'mix_file': '2014_08_05 Rezeptur_MI.xlsx'}

    with open(metadataPath, 'w') as yamlFile:
        yaml.dump(target_data, yamlFile)
    print("Generate knowledgeGraph:")
    # generate_knowledge_graph(metadataPath, filename)

    # g = rdflib.Graph()
    # result = g.parse(filename, format='ttl')
    # print("Query knowledgeGraph:")
    #
    # print_first()
    #
    # result = query_and_test_result('diameter', 'ns1:informationbearingentity1',
    #                                'ns7:has_decimal_value', g)
    # for row in result:
    #     assert (target_data['diameter'] = float(row.o))
    #
    # print_next()
    #
    # result = query_and_test_result('length', 'ns1:informationbearingentity2',
    #                                'ns7:has_decimal_value', g)
    # for row in result:
    #     assert(target_data['length'] = float(row.o))
    #
    # print_next()
    #
    # result = query_and_test_result('weight', 'ns1:informationbearingentity3',
    #                                'ns7:has_decimal_value', g)
    # for row in result:
    #     assert(target_data['weight'] = float(row.o))
    #
    # print_next()
    #
    # time = get_date_time_value(target_data)
    # result = query_and_test_result('operator_time', 'ns1:informationbearingentity6',
    #                                'ns1:has_datetime_value', g)
    # for row in result:
    #     assert(time, datetime.datetime.strptime(str(row.o).replace('T', ' '), "%Y-%m-%d %H:%M:%S"))
    #
    # print_next()
    #
    # result = query_and_test_result('operatorName', 'ns1:informationbearingentity7',
    #                                'ns7:has_text_value', g)
    # for row in result:
    #     assertEqual(target_data['tester_name'] = str(row.o))
    #
    # print_next()
    # result = query_and_test_result('rawFileName', 'ns1:informationbearingentity4',
    #                                'ns7:has_text_value', g)
    # for row in result:
    #     assert(target_data['experimentName'] = str(row.o))
    #
    # print_next()
    # result = query_and_test_result('rawFileName', 'ns1:informationbearingentity5',
    #                                'ns7:has_text_value', g)
    # rawPath = os.path.join(Path(__file__).parents[2], 'usecases',
    #                        'MinimumWorkingExample', 'Data', 'E-modul', target_data['experimentName'])
    #
    # for row in result:
    #     assert(rawPath = str(row.o))
    #
    # print('Pass all tests')
    #
    # print("Deleting test files:")
    # os.remove(filename)
    # os.remove(metadataPath)
    #
    # print("End testing")

# test_generate_knowledge_graph()
