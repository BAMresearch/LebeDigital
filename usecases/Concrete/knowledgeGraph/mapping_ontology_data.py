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


# In[2]:


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
COM = My_world.get_namespace('https://github.com/BAMresearch/ModelCalibration/blob/Datasets/usecases/Concrete/ConcreteOntology/Concrete_Ontology_MSEO.owl')
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
    float(data['remark'][i].split()[1].replace(',','.')) for i in data.index
]
data['control unit'] = [
    data['remark'][i].split()[2] for i in data.index
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
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_')).iri)), 
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
        # add testers (prüfers) as instances in class Agent
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Agent(data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Agent.iri))
            )
        )
        # add person/tester in class DesignativeName
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.DesignativeName(data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.DesignativeName.iri))
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
        # add unit in MeasurementUnitOfForceRate
        g.add(
            (
                URIRef(urllib.parse.quote(COM.MeasurementUnitOfForceRate(data['control unit'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(COM.MeasurementUnitOfForceRate.iri))
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
                URIRef(urllib.parse.quote(CCO.Mass(data['weight'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.Mass.iri))
            )
        )
#        
        # add specimen region as instances in class MeasurementRegion
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(lebedigital_concrete.MeasurementRegion.iri))
            )
        )
        # add weight, diameter, length, dataset path, tester name value, force rate value
        # in class InformationBearingEntity
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['weight'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['diameter'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['length'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['sample name'][i].replace(' ','_') + '_' + data['control'][i] ).iri)), 
                RDF.type, 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity.iri))
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

# In[20]:


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
        # Diameter, Length, Mass, ForceRate, RawDataSet, DesignativeName obo:RO_0010001 InformationBearingEntity
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Diameter(data['diameter'][i]).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['diameter'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Diameter(data['length'][i]).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['length'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.Diameter(data['weight'][i]).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['weight'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(ConcreteMSEO_ontology.ForceRate(data['control'][i]).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['control'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(lebedigital_concrete.RawDataSet(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.DesignativeName(data['tester'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(OBO.RO_0010001.iri)),
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['tester'][i].replace(' ','_')).iri))
            )
        )
        # InformationBearingEntity of ForceRate uses_measurement_unit MeasurementUnitOfForceRate
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['sample name'][i].replace(' ','_') + '_' + data['control'][i] ).iri)), 
                URIRef(urllib.parse.quote(CCO.uses_measurement_unit.iri)), 
                URIRef(urllib.parse.quote(COM.MeasurementUnitOfForceRate(data['control unit'][i]).iri))
            )
        )


# <h5 style="color:#1f5dbf">adding data with data properties</h5>  

# In[21]:


# add weight, diameter, length, dataset path, tester name value, force rate value
for i in data.index:
    with lebedigital_concrete:
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['weight'][i]).iri)), 
                URIRef(urllib.parse.quote(CCO.has_decimal_value.iri)),
                Literal(data['weight_number'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['diameter'][i]).iri)), 
                URIRef(urllib.parse.quote(CCO.has_decimal_value.iri)),
                Literal(data['diameter_number'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['length'][i]).iri)), 
                URIRef(urllib.parse.quote(CCO.has_decimal_value.iri)),
                Literal(data['length_number'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['sample name'][i].replace(' ','_') + 'specimen.dat').iri)), 
                URIRef(urllib.parse.quote(CCO.has_URI_value.iri)), 
                Literal(data['file path'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['tester'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.quote(CCO.has_text_value.iri)), 
                Literal(data['tester'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(CCO.InformationBearingEntity(data['sample name'][i].replace(' ','_') + '_' + data['control'][i]).iri)), 
                URIRef(urllib.parse.quote(CCO.has_decimal_value.iri)), 
                Literal(data['control value'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.quote(COM.MeasurementUnitOfForceRate(data['control unit'][i]).iri)), 
                URIRef(urllib.parse.quote(CCO.has_text_value.iri)), 
                Literal(data['control unit'][i])
            )
        )


g.serialize(destination=graphPath, format="turtle")


# In[30]:

print('---------------------------------------------------------------------')
print('sample queries')
print('number of EModul experiments: ', 39)
print('number of classes in Emodul ontology: ', len(list(lebedigital_concrete.classes())))

q = """
    

    SELECT (count(*) as ?Triples)
    WHERE {
        {?s ?p ?o}
    }
"""


# In[31]:


results = g.query(q)
for result in results:
    print('number of triples', result)


# In[ ]:




