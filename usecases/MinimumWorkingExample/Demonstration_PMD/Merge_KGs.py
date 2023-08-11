from rdflib import Graph, URIRef

def merge_graphs(p_graph1, p_graph2, p_graph3):
    merged_graph = Graph()

    graph1 = Graph()
    graph1.parse(p_graph1,
                 format="turtle")  # Load RDF data from a file

    graph2 = Graph()
    graph2.parse(p_graph2,
                 format="turtle")  # Load another RDF data from a file

    graph3 = Graph()
    graph3.parse(p_graph3,
                 format="turtle")  # Load another RDF data from a file

    for s, p, o in graph1:
        merged_graph.add((s, p, o))

    for s, p, o in graph2:
        merged_graph.add((s, p, o))

    for s, p, o in graph3:
        merged_graph.add((s, p, o))

    fileName = "MergedGraph.ttl"
    with open(fileName, "w") as merged_file:
        merged_file.write(merged_graph.serialize(format="turtle"))

    return fileName

def main():
    # Create two example RDF graphs
    #graph1 = Graph()
    #graph1.parse("../../../usecases/MinimumWorkingExample/Mapping_Example/testMixMapped.ttl", format="turtle")  # Load RDF data from a file

    #graph2 = Graph()
    #graph2.parse("../../../usecases/MinimumWorkingExample/Mapping_Example/testSpecimenMapped.ttl", format="turtle")  # Load another RDF data from a file

    graph1 = "../../../usecases/MinimumWorkingExample/Mapping_Example/testMixMapped.ttl"
    graph2 = "../../../usecases/MinimumWorkingExample/Mapping_Example/testSpecimenMapped.ttl"
    graph3 = "../../../usecases/MinimumWorkingExample/Mapping_Example/testMapped.ttl"
    # Merge the RDF graphs using the merge_graphs function
    merged_graph = merge_graphs(graph1, graph2, graph3)


    # Print merged graph triples
    #for s, p, o in merged_graph:
     #   print(s, p, o)

if __name__ == "__main__":
    main()
