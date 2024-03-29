# standard library
import os
import yaml
import datetime

# third party imports
import numpy as np
import arviz as az
from owlready2 import *
from yaml.loader import SafeLoader
from rdflib import URIRef, Graph, Literal, BNode
from rdflib.namespace import FOAF, RDF, XSD
from pathlib import Path

#This script is fully dependent on the ontology that on can find at
#lebedigital/ConcreteOntology/ 
#Here we generate the Path to get all the ontologies needed to run the script
#Warnings: If the architecture of the lebedigital change this may need to be changed

#Settings for files and directories:
ontologyDirName = 'ConcreteOntology'
lebedigitalPkgName = 'lebedigital'
baseDir = Path(__file__).parents[2]
ontologyPath = os.path.join(baseDir, lebedigitalPkgName, ontologyDirName)
MinWorkingExpDir = os.path.join(baseDir, 'usecases', 'MinimumWorkingExample')
rawBaseDir = os.path.join(MinWorkingExpDir, 'Data', 'E-modul')
processedDataDir =os.path.join(MinWorkingExpDir, 'emodul', 'processed_data')


def import_metadata(locationOfMetadata):
    """Import metadata from a file and return it as a dict
    
    Parameters
    ----------
    locationOfMetadata: string
        Path to the metadata File
    
    Returns
    -------
    metaDataDict: dict
        Dictionary that represent the metadata in the given file"""
    with open(locationOfMetadata) as metaDataFile:
        try:
            metaDataDict = yaml.load(metaDataFile, Loader=SafeLoader)
            return metaDataDict
        except Exception as e:
            print("Path error: ", file=sys.stderr)
            print(e, file=sys.stderr)

###############################################################################
#Utilitys that can be useful for testing if we need to change the script
def print_ontology_classes(onto):
    c_list = list(onto.classes())
    print("List of classes: ")
    for c in c_list:
        print(c)
    print(f"There are '{len(c_list)}' classes")

def print_ontology_instances(onto):
    i_list = list(onto.individuals())
    print("List of individuals: ")
    for i in i_list:
        print(i)
    print(f"There are '{len(i_list)}' instancies")

###############################################################################
def export_knowledge_graph_to_ttl(onto, locationOfKnowledgeGraph):
    """Export the knowledge graph that we generate by adding instances to our ontology
    Parameters
    ----------
    onto: owlready2 ontology
        the object containing our knowledge graph
    locationOfKnowledgeGraph: str
        Path to the file
    
    Returns:
    --------
        void

    Warnings:
    -Since owlready2 is not used to create knowledge graphs 
    we need to export it using RDFLIB.
    -Owlready2 is not designed for knowledge graph generation so using it this 
    way may result in ambiguous results.

    """
    knowledge_graph = onto.as_rdflib_graph()
    knowledge_graph.serialize(destination=locationOfKnowledgeGraph, format="turtle")

def make_valid_uri_from_string(s):
    return s.replace('.','_').replace(' ','_').replace(':','_').replace(',','_')

def get_date_time_value(metadata):
    str_date_time = metadata['operator_date'] + ' ' +metadata['operator_timestamp']
    date_time = datetime.datetime.strptime(str_date_time, "%d.%m.%Y %H:%M:%S")
    return date_time
    
def generate_knowledge_graph(metadataPath, knowledgeFile):
    """
    Generate a knowledge graph

    Parameters
    ----------
    metadataPath: str
        Path to the file containing the metadata
    knowledgeFile: str
        Path to the file in which we want to export the knowledge graph
    
    Returns
    -------
        void
    """    
    My_world = World()

    #Load all ontologies
    cco_ontology = My_world.get_ontology(ontologyPath + "/MergedAllCoreOntology.owl").load()
    mseo_mid = My_world.get_ontology(ontologyPath + "/MSEO_mid.owl").load()
    PeriodicTable_ontology = My_world.get_ontology(ontologyPath + "/PeriodicTable.owl").load()
    WCTmidonto = My_world.get_ontology(ontologyPath + "/WCTmid.owl").load()
    CSTonto = My_world.get_ontology(ontologyPath + "/ConcreteStressTestOntologie.owl").load()
    impO = My_world.get_ontology(os.path.join(ontologyPath,'EM.xml')).load()
    ConcreteMSEO_ontology = My_world.get_ontology(ontologyPath + "/Concrete_Ontology_MSEO.owl").load()

    #Get all namespaces
    bfo = impO.get_namespace("http://purl.obolibrary.org/obo/")
    cco = impO.get_namespace("http://www.ontologyrepository.com/CommonCoreOntologies/")
    mseo = impO.get_namespace("https://purl.matolab.org/mseo/mid/")
    con = impO.get_namespace('http://w3id.org/concrete/#')
    xml = impO.get_namespace('http://www.w3.org/2001/XMLSchema#')

    #I need to redefined the property because the property was not defined as functional
    with impO:
        # infoBearing = cco.InformationBearingEntity
        # class has_decimal_value(DataProperty):
        #     domain    = [cco.InformationBearingEntity]
        #     range     = [float]
        class has_text_value(DataProperty, FunctionalProperty):
            domain    = [cco.InformationBearingEntity]
            range     = [str]
    #Extract metadata
    with open(metadataPath) as metadataFile:
        try:
            metadata = yaml.load(metadataFile, Loader=SafeLoader)
            rawPath = os.path.join(rawBaseDir, metadata['experimentName'])
            processedFile = os.path.join(processedDataDir, metadata['experimentName'])
            ###########################################################################
            ########################ADD INDIVIDUALS####################################
            ###########################################################################
            
            #Root
            concreteSpecimenIndividual = impO.Specimen(make_valid_uri_from_string(metadata['experimentName']))

            #Individuals needed to give a Diameter value
            specimenDiameter = cco.Diameter()
            specimenDiameterValue = cco.InformationBearingEntity()
            
            #Individuals for length
            specimenLength = cco.Length()
            specimenLengthValue = cco.InformationBearingEntity()
            
            #Individuals for Mass
            specimenMass = cco.Mass()
            specimenMassValue = cco.InformationBearingEntity()

            #Add comment'
            specimenSecantModulus = ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(make_valid_uri_from_string(metadata['experimentName']))

            #Add filename, filePath of RawData
            specimenRawDatafile = mseo.RawDataSet()
            specimenRawFilename =cco.InformationBearingEntity()
            specimenRawFilePath = cco.InformationBearingEntity()

            #Add date
            experimentDate = cco.Day()
            experimentDateValue = cco.InformationBearingEntity()

            #Add operator Name
            experimentOperator = cco.Agent()
            experimentOperatorName = cco.DesignativeName()
            experimentOperatorNameValue = cco.InformationBearingEntity()

            #Add Filename, filePath of Processed Data
            specimenProcessingOfData = bfo.BFO_0000015()
            specimenProcessedDatafile = mseo.AnalysedDataSet()
            specimenProcessedFilename = cco.InformationBearingEntity()
            specimenProcessedPath = cco.InformationBearingEntity()
            
            ###########################################################################
            #########################ADD PROPERTIES####################################
            ###########################################################################
            
            #Begin Path to SecantDetermination
            concreteSpecimenIndividual.is_input_of= [specimenSecantModulus]

            #Add path to diameter value
            concreteSpecimenIndividual.RO_0000086 = [specimenDiameter]
            specimenDiameter.RO_0010001 = [specimenDiameterValue]
            specimenDiameterValue.has_decimal_value = [float(metadata['diameter'])]
            specimenDiameterValue.uses_measurement_unit = [cco.MillimeterMeasurementUnit]
        
            #Add path to length value
            concreteSpecimenIndividual.RO_0000086 =[specimenLength]
            specimenLength.RO_0010001 = [specimenLengthValue]
            specimenLengthValue.has_decimal_value = [float(metadata['length'])]
            #Add length unit
            if(metadata['length_unit'] == "mm"):
                specimenLengthValue.uses_measurement_unit = [cco.MillimeterMeasurementUnit]
            
            #Add path to Mass
            concreteSpecimenIndividual.RO_0000086 = [specimenMass]
            specimenMass.RO_0010001 = [specimenMassValue]
            specimenMassValue.has_decimal_value = [float(metadata['weight'])]
            #Add unit
            if(metadata['weight_unit'] == 'g'):
                specimenMassValue.uses_measurement_unit = [cco.GramMeasurementUnit]

            #Add path to Rawfile
            specimenSecantModulus.has_output = [specimenRawDatafile]
            specimenRawDatafile.RO_0010001 = [specimenRawFilename, specimenRawFilePath]
            specimenRawFilename.has_text_value =metadata['experimentName']
            specimenRawFilePath.has_text_value = rawPath

            #Add path to Processedfile
            specimenRawDatafile.is_input_of = [specimenProcessingOfData]
            specimenProcessingOfData.has_output = [specimenProcessedDatafile]
            specimenProcessedDatafile.RO_0010001 = [specimenProcessedFilename,
                                                    specimenProcessedPath]
            specimenProcessedFilename.has_text_value = metadata['experimentName'] + '.csv'
            specimenProcessedPath.has_text_value = processedDataDir
            #Add path to date
            specimenSecantModulus.occures_on = [experimentDate]
            experimentDate.RO_0010001 = [experimentDateValue]
            operatorTime = get_date_time_value(metadata)
            experimentDateValue.has_datetime_value = [operatorTime]

            #Add Path to operator
            specimenSecantModulus.has_agent = [experimentOperator]
            experimentOperator.designated_by = [experimentOperatorName]
            experimentOperatorName.RO_0010001 =[experimentOperatorNameValue]
            experimentOperatorNameValue.has_text_value = metadata['tester_name']

            ###########################################################################
            #########################Export Graph######################################
            ###########################################################################
            
            export_knowledge_graph_to_ttl(My_world, knowledgeFile)
        except Exception as e:
            print("Path error: ", file=sys.stderr)
            print(e, file=sys.stderr)
    

# metadataPath = "/home/gilif/BAM/LeBeDigital_Projects/mapping_script_Lebedigital/usecases/MinimumWorkingExample/emodul/metadata_yaml_files/testMetaData.yaml"

# generate_knowledge_graph(owlPath, metadataPath)