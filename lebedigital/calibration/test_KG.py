import rdflib
def query_objects(queries, graph):
        # function to get objects from specific subject, predicate pairs, knowing there is only one result
        results = {}
        for query in queries:
                q = f"""
                        SELECT ?object
                        WHERE {{
                                {queries[query]['subject']}  {queries[query]['predicate']} ?object
                        }}
                """

                result = graph.query(q)
                for r in result:
                        results[query] = r["object"]

        return results

# define queries
queries = {
        'diameter' :  {'subject' : 'ns1:informationbearingentity1', 'predicate' : 'ns1:has_decimal_value' },
        'length' :    {'subject' : 'ns1:informationbearingentity2', 'predicate' : 'ns1:has_decimal_value' },
        'path' :      {'subject' : 'ns1:informationbearingentity9', 'predicate' : 'ns9:has_text_value' },
        'file_name' : {'subject' : 'ns1:informationbearingentity8', 'predicate' : 'ns9:has_text_value' }
}

# Path to KG
path_to_KG = '../../usecases/MinimumWorkingExample/emodul/knowledge_graphs/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.tll'

# initialize the graph
knowledge_graph = rdflib.Graph()
knowledge_graph.parse(path_to_KG, format='ttl')

# queries
results = query_objects(queries,knowledge_graph)

# geometry
length = results['length']
print('Length: ', length)
diameter = results['diameter']
print('Diameter: ',diameter)

# experimental data
file_path = results['path'] + '/' + results['file_name']
print(file_path)