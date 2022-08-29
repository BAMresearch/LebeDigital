import pytest
import rdflib
import datetime
import os
import yaml

from pathlib import Path

from lebedigital.mapping.emodul_mapping import generate_knowledge_graph, get_date_time_value

def calibrate_E_from_KG(KG_file):
    print('Path to file', KG_file)

    g = rdflib.Graph()
    result = g.parse(KG_file, format='ttl')
    print("Query knowledgeGraph:")

    result = query_and_test_result('diameter', 'ns1:informationbearingentity1','ns7:has_decimal_value', g)
    print(result)
    for row in result:
        print(row)
        print('TEST')
        print('Diameter:',float(row.o))


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

path_to_KG = '../../usecases/MinimumWorkingExample/emodul/knowledge_graphs/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.tll'
calibrate_E_from_KG(path_to_KG)