import pytest
import rdflib
import datetime
import os
import yaml

from pathlib import Path

#Path to KG
path_to_KG = '../../usecases/MinimumWorkingExample/emodul/knowledge_graphs/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.tll'

#This Create the graph
g = rdflib.Graph()

#Now we parse the file and add to the graph
#So now g contains the graph
g.parse(path_to_KG, format='ttl') #We need to declare the format because the default is rdf

#Example of query
#This query allows to get the diameter
#In natural language it says:
#Give me all objects that are decimal value of diameter
q = """
        SELECT ?object
        WHERE {
                ns1:informationbearingentity1 ns1:has_decimal_value ?object
        }
"""
#I know that ns1:informationbearingentity1 isn't meaningful but it's the name
#given by owlready2 when I declare a diameter. I can change the name but before 
#doing that we need to decide on a standard.
#The name is derived from the ontology for the moment. If they get better labels
# We will be able to have more meaningful names.

#You can quickly find the name the entity you want by looking 
#Query and store the result
result = g.query(q)

#g.query return an iterable containing objects that representes each result.
#Thats why I use a loop to read the result.
#The documentation say that it's the best way to go with them
for r in result:
        print("Diameter: " + r["object"])

#The keys are the variables after ?
#More examples:
print("############New Results#############")
#Get all entities that have 98.6 as diameter
q = """
        SELECT ?predicate
        WHERE {
                ?predicate ns1:has_decimal_value 98.6
        }
"""

result = g.query(q)
for r in result:
        print("Diameter: " + r["predicate"])

q = """
        SELECT ?predicate ?prop
        WHERE {
                ?predicate ?prop 98.6
        }
"""
print("############New Results#############")
result = g.query(q)
for r in result:
        print("Predicate: " + r["predicate"])
        print("Prop: " + r["prop"])

#Get the weigth:

q = """
        SELECT ?obj
        WHERE {
                ns1:informationbearingentity2 ns1:has_decimal_value ?obj
        }
"""
print("############New Results#############")
result = g.query(q)
for r in result:
        print("Weigth: " + r["obj"])