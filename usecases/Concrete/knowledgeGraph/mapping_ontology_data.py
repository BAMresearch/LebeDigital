#!/usr/bin/env python
# coding: utf-8

# 
# <h1 style="color:#fc6b03">Mapping metadata and ontology</h1>
# 
# ### Data and metadata
# - Data from <span style="color:orange">Young modulus</span> experiment of concrete material in .dat file (table format)
# - Metadata was extracted from the header of the file (listed all metadata in a csv file)
# 
# ### Ontology
# - The <span style="color:orange">ontology BWMD_mid</span> as mid ontology
# - <span style="color:orange">Concrete_Ontology_rdf</span> as specific usecase in LeBeDigital
# 
# ### The workflow of the script
# 1. Installing and importing the necessary libraries: <span style="color:orange">owlready, rdflib and pandas</span>
# 2. Using owlready to load the 2 ontology files and pandas to load csv file
# 3. Using rdflib to create a graph and add triples by adding data from csv file and ontology entities from ontology files
# 4. Exporting the graph as .ttl file (turtle file)
# 

# In[1]:


from owlready2 import *
import pandas as pd
import urllib.parse
from pathlib import Path
import os

# In[23]:

baseDir = Path(__file__).resolve().parents[1]
ontologyPath = os.path.join(baseDir,'ConcreteOntology')
metadataPath = 'data/metadata.csv'
graphPath = 'data/EM_Graph.ttl'
importedOntologiesPath = os.path.join(baseDir,'GraphCreation/Ontologies')


# <h3 style="color:#1f5dbf">load concrete material ontology for Emodul experiment</h3>   

# In[3]:


onto_path.append(".")
My_world = World()

cco_ontology = My_world.get_ontology(importedOntologiesPath + "/MergedAllCoreOntology.owl").load()
mseo_mid = My_world.get_ontology(importedOntologiesPath + "/MSEO_mid.owl").load()
PeriodicTable_ontology = My_world.get_ontology(importedOntologiesPath + "/PeriodicTable.owl").load()
WCTmidonto = My_world.get_ontology(importedOntologiesPath + "/WCTmid.owl").load()
CSTonto = My_world.get_ontology(importedOntologiesPath + "/ConcreteStressTestOntologie.owl").load()
lebedigital_concrete = My_world.get_ontology(os.path.join(ontologyPath,'EM.xml')).load()
ConcreteMSEO_ontology = My_world.get_ontology(ontologyPath + "/Concrete_Ontology_MSEO.owl").load()

BWMD = My_world.get_namespace("https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#")
WCTmid = My_world.get_namespace("https://mobi.com/ontologies/6/2021/WCTmid#")
#WCT = My_world.get_namespace("https://mobi.com/ontologies/6/202https://mobi.com/ontologies/6/2021/WoodCompressionTest#1/WoodCompressionTest#")
CCO = My_world.get_namespace("http://www.ontologyrepository.com/CommonCoreOntologies/")
CST = My_world.get_namespace("https://mobi.com/ontologies/7/2021/ConcreteStressTestOntologie#")
COM = My_world.get_namespace('https://git.material-digital.de/lebedigital/concreteontology/-/blob/Ontology/Concrete_Ontology_MSEO.owl')
OBO = My_world.get_namespace("http://purl.obolibrary.org/obo/")

lebedigital_concrete.imported_ontologies.append(cco_ontology)
lebedigital_concrete.imported_ontologies.append(PeriodicTable_ontology)
lebedigital_concrete.imported_ontologies.append(WCTmidonto)
lebedigital_concrete.imported_ontologies.append(CSTonto)
lebedigital_concrete.imported_ontologies.append(mseo_mid)


# <h3 style="color:#1f5dbf">metadata from Emodul experiments</h3>  

# In[7]:


data = pd.read_csv(metadataPath)


# In[8]:


data.head()


# In[9]:


data['experiment name'] = ['E-modul experiment ' + data['sample name'][i]
                           for i in data.index
                          ]


# In[10]:


data.head()


# <h5 style="color:#1f5dbf">splitting string from remark in control and the value of control force/stress</h5>  

# In[11]:


data['control'] = [
    data['remark'][i].split()[0] for i in data.index
]
data['control value'] = [
    data['remark'][i].split()[1] + data['remark'][i].split()[2] for i in data.index
]


# In[12]:


data.head()


# <h5 style="color:#1f5dbf">add 1 column for example file path in data</h5>  

# In[13]:


data['file path'] = [
    'https://github.com/BAMresearch/ModelCalibration/tree/main/usecases/Concrete/Data/E-modul'
    + '/' + data['sample name'][i]
    for i in data.index
]


# In[14]:


data.head()


# <h5 style="color:#1f5dbf">convert string to number in the columns</h5>  

# In[15]:


data['weight_number'] = [
    float(data['weight'][i].replace(',','.'))
    for i in data.index
]
data['diameter_number'] = [
    float(data['diameter'][i].replace(',','.'))
    for i in data.index
]
data['length_number'] = [
    float(data['length'][i].replace(',','.'))
    for i in data.index
]


# In[16]:


data.head()


# <h3 style="color:#1f5dbf">Generate a graph and add triples to this graph</h3>  

# In[17]:


from rdflib import URIRef, Graph, Literal, BNode
from rdflib.namespace import FOAF, RDF
g = My_world.as_rdflib_graph()


# <h5 style="color:#1f5dbf">adding instances into classes</h5>  

# In[18]:


for i in data.index:
    with lebedigital_concrete:
        # add experiments as instances in class DeterminationOfSecantModulusOfElasticity
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity.iri))
            )
        )
        # add specimens as instances in class Specimen
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen(data['sample name'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen.iri))
            )
        )
        # add testers (prüfers) as instances in class Person
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Agent(data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Agent.iri))
            )
        )
        # add start date as instances in class Day
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Day(data['operator date'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Day.iri))
            )
        )
        # add controls as instances in class ForceRate
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.ForceRate(data['control'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.ForceRate.iri))
            )
        )
        # add length of the specimen in class Length
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Length(data['length'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Length.iri))
            )
        )
        # add Diameter of the specimen in class Diameter
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Diameter(data['diameter'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Diameter.iri))
            )
        )
        # add Weight of the specimen in class Mass
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Mass(data['diameter'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Mass.iri))
            )
        )
#         # add name of the operator in class DesignativeName
#         g.add(
#             (
#                 Literal(data['tester'][i]), 
#                 RDF.type, 
#                 URIRef(urllib.parse.quote(CCO.DesignativeName.iri))
#             )
#         )
        # add specimen region as instances in class MeasurementRegion
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(lebedigital_concrete.MeasurementRegion.iri))
            )
        )


# <h5 style="color:#1f5dbf">adding data type in process dataset</h5>  

# In[19]:


for i in data.index:
    with lebedigital_concrete:
        g.add(
                (
                    URIRef(urllib.parse.quote(lebedigital_concrete.RawDataSet(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.quote(lebedigital_concrete.RawDataSet.iri))
                )
            )


# <h5 style="color:#1f5dbf">adding data with object properties and data properties</h5>  

# In[21]:


for i in data.index:
    with lebedigital_concrete:
        # experiment has_output RawDataSet ('specimen.dat')
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri)), 
                URIRef(urllib.parse.quote(CCO.has_output.iri)), 
                URIRef(urllib.parse.quote(lebedigital_concrete.RawDataSet(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri))
            )
        )
        # experiment hasOperator from class person
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri)), 
                URIRef(urllib.parse.quote(CCO.has_agent.iri)), 
                URIRef(urllib.parse.quote(CCO.Agent(data['tester'][i].replace(' ','_')).iri))
            )
        )
        # specimen is_input_of experiment
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen(data['sample name'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(CCO.is_input_of.iri)), 
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri))
            )
        )
        # ForceRate is_input_of experiment
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.ForceRate(data['control'][i]).iri)), 
                URIRef(urllib.parse.quote(CCO.is_input_of.iri)), 
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri))
            )
        )
        # Experiment occures_on Day
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri)), 
                URIRef(urllib.parse.quote(CCO.occures_on.iri)), 
                URIRef(urllib.parse.quote(CCO.Day(data['operator date'][i]).iri))
            )
        )
        # specimen obo:BFO_0000051 MeasurementRegion
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen(data['sample name'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(OBO.BFO_0000051.iri)), 
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name'][i].replace(' ','_')).iri))
            )
        )
        # MeasurementRegion obo:BFO_0000086 Diameter, Length, Mass
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0000086.iri)), 
                URIRef(urllib.parse.quote(CCO.Diameter(data['diameter'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0000086.iri)), 
                URIRef(urllib.parse.quote(CCO.Diameter(data['length'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0000086.iri)), 
                URIRef(urllib.parse.quote(CCO.Diameter(data['weight'][i]).iri))
            )
        )


# <h5 style="color:#1f5dbf">create more data properties </h5>  

# with bwmd_mid:
#     class hasSpecimenName(DataProperty):
#         range = [str]
#     class hasName(DataProperty):
#         range = [str]
#     class hasDiameter(DataProperty):
#         range = [str]
#     class hasWeight(DataProperty):
#         range = [str]
#     class hasLength(DataProperty):
#         range = [str]
# with lebedigital_concrete:
#     class hasValue(DataProperty):
#         range = [str]
#     class hasFilePath(DataProperty):
#         range = [str]
#     

# <h5 style="color:#1f5dbf">data properties for specimen and remark of experiment </h5>  

# for i in data.index:
#     # specimen hasSpecimenName
#     g.add(
#         (
#             URIRef(urllib.parse.quote(bwmd_mid.BWMD_00048(data['sample name'][i].replace(' ','_')).iri)), 
#             URIRef(urllib.parse.quote(bwmd_mid.hasSpecimenName.iri)),
#             Literal(data['sample name'][i].replace(' ','_'))
#         )
#     )
#     # specimen hasDiameter
#     g.add(
#         (
#             URIRef(urllib.parse.quote(bwmd_mid.BWMD_00048(data['sample name'][i].replace(' ','_')).iri)), 
#             URIRef(urllib.parse.quote(bwmd_mid.hasDiameter.iri)),
#             Literal(data['diameter'][i])
#         )
#     )
#     # specimen hasWeight
#     g.add(
#         (
#             URIRef(urllib.parse.quote(bwmd_mid.BWMD_00048(data['sample name'][i].replace(' ','_')).iri)), 
#             URIRef(urllib.parse.quote(bwmd_mid.hasWeight.iri)),
#             Literal(data['weight'][i])
#         )
#     )
#     # specimen hasLength
#     g.add(
#         (
#             URIRef(urllib.parse.quote(bwmd_mid.BWMD_00048(data['sample name'][i].replace(' ','_')).iri)), 
#             URIRef(urllib.parse.quote(bwmd_mid.hasLength.iri)),
#             Literal(data['length'][i])
#         )
#     )
#     # StressRate hasValue
#     g.add(
#         (
#             URIRef(urllib.parse.quote(lebedigital_concrete.StressRate(data['control'][i]).iri)), 
#             URIRef(urllib.parse.quote(lebedigital_concrete.hasValue.iri)), 
#             Literal(data['control value'][i])
#         )
#     )
#     # file hasFilePath
#     g.add(
#         (
#             URIRef(urllib.parse.quote(bwmd_mid.BWMD_00068(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri)), 
#             URIRef(urllib.parse.quote(lebedigital_concrete.hasFilePath.iri)), 
#             Literal(data['file path'][i])
#         )
#     )

# <h5 style="color:#1f5dbf">data properties for person </h5>  

# for i in data.index:
#     # person hasName
#     g.add(
#         (
#             URIRef(urllib.parse.quote(bwmd_mid.BWMD_00004(data['tester'][i].replace(' ','_')).iri)), 
#             URIRef(urllib.parse.quote(bwmd_mid.hasName.iri)), 
#             Literal(data['tester'][i].replace(' ','_'))
#         )
#     )

# In[24]:


g.serialize(destination=graphPath, format="turtle")


# In[26]:


q = """
    prefix ns1: <http://www.w3.org/2002/07/owl#> 
    prefix ns10: <https://www.materials.fraunhofer.de/ontologies/graph_designer#> 
    prefix ns11: <http://purl.org/dc/terms/> 
    prefix ns2: <http://www.ontologyrepository.com/CommonCoreOntologies/> 
    prefix ns3: <http://purl.obolibrary.org/obo/> 
    prefix ns4: <http://www.geneontology.org/formats/oboInOwl#> 
    prefix ns5: <http://www.daml.org/2003/01/periodictable/PeriodicTable#> 
    prefix ns6: <http%3A//www.ontologyrepository.com/CommonCoreOntologies/> 
    prefix ns7: <http%3A//purl.obolibrary.org/obo/> 
    prefix ns8: <http://purl.org/dc/elements/1.1/> 
    prefix ns9: <https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#> 
    prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
    prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    prefix xsd: <http://www.w3.org/2001/XMLSchema#> 
    prefix con: <https://git.material-digital.de/lebedigital/concreteontology/-/blob/Concrete_Ontology_MSEO.owl>

    SELECT ?o
    WHERE {
        {
        <https%3A//purl.matolab.org/mseo/mid/E-modul_experiment_BA-Losert_E-Modul_28d_v._04.08.14_Probe_4>
        ns6:has_agent 
        ?o
        }
        UNION
        {
        <https%3A//purl.matolab.org/mseo/mid/E-modul_experiment_BA-Losert_E-Modul_28d_v._04.08.14_Probe_4>
        ns6:occures_on 
        ?o
        }
    }
"""


# In[28]:


results = g.query(q)


# In[31]:

print('Example query results for the question: when the experiment E-modul_experiment_BA-Losert_E-Modul_28d_v._04.08.14_Probe_4  occured and by whom?')
for result in results:
    print(result)


# In[ ]:




