from owlready2 import *
import yaml
from yaml.loader import SafeLoader
import urllib.parse
from pathlib import Path
import os
import datetime
from rdflib import URIRef, Graph, Literal, BNode
from rdflib.namespace import FOAF, RDF


BASEDIR4 = Path(__file__).resolve().parents[4]
BASEDIR1 = Path(__file__).resolve().parents[1]
BASEDIR2 = Path(__file__).resolve().parents[2]
E_MODUL_RAWDATA_PATH = os.path.join(os.path.join(os.path.join(BASEDIR2,'Example'),'Data'),'E-modul')
E_MODUL_PROCESSED_DATA_PATH = os.path.join(os.path.join(os.path.join(BASEDIR2,'Example'),'emodul'),'processeddata')
ontologyPath = os.path.join(BASEDIR4, 'lebedigital', 'ConcreteOntology')

onto_path.append(".")
My_world = World()

cco_ontology = My_world.get_ontology(ontologyPath + "/MergedAllCoreOntology.owl").load()
mseo_mid = My_world.get_ontology(ontologyPath + "/MSEO_mid.owl").load()
PeriodicTable_ontology = My_world.get_ontology(ontologyPath + "/PeriodicTable.owl").load()
WCTmidonto = My_world.get_ontology(ontologyPath + "/WCTmid.owl").load()
CSTonto = My_world.get_ontology(ontologyPath + "/ConcreteStressTestOntologie.owl").load()
lebedigital_concrete = My_world.get_ontology(os.path.join(ontologyPath,'EM.xml')).load()
ConcreteMSEO_ontology = My_world.get_ontology(ontologyPath + "/Concrete_Ontology_MSEO.owl").load()

BWMD = My_world.get_namespace("https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#")
WCTmid = My_world.get_namespace("https://mobi.com/ontologies/6/2021/WCTmid#")
CCO = My_world.get_namespace("http://www.ontologyrepository.com/CommonCoreOntologies/")
CST = My_world.get_namespace("https://mobi.com/ontologies/7/2021/ConcreteStressTestOntologie#")
CON = My_world.get_namespace('http://w3id.org/concrete')
OBO = My_world.get_namespace("http://purl.obolibrary.org/obo/")
My_world.get_namespace('https://purl.matolab.org/mseo/mid/')

lebedigital_concrete.imported_ontologies.append(cco_ontology)
lebedigital_concrete.imported_ontologies.append(PeriodicTable_ontology)
lebedigital_concrete.imported_ontologies.append(WCTmidonto)
lebedigital_concrete.imported_ontologies.append(CSTonto)
lebedigital_concrete.imported_ontologies.append(mseo_mid)

def make_valid_uri_from_string(s):
    return s.replace('.','_').replace(' ','_').replace(':','_').replace(',','_')

def string_to_number(s):
    if type(s) == int or type(s) == float:
        return s
    else:
        return float(s.replace(',','.'))

def metadata_ontology_mapping(locationOfMetadata, locationOfKnowledgeGraph):
    try:
        with open(locationOfMetadata) as f:
            data = yaml.load(f, Loader=SafeLoader)
        nameOfExperiment = data['experimentName']
        operatorTime = datetime.datetime.strptime(data['operatorTimestamp'], '%d.%m.%Y %H:%M:%S')
        operatorDay = data['operatorDate']
        operator = make_valid_uri_from_string(data['tester'])
        dataType = data['dataType']
        specimenName = data['specimenName']
        # convert string to number of some specimen parameters values
        specimenDiameter = string_to_number(data['diameter'])
        specimenLength = string_to_number(data['length'])
        specimenMass = string_to_number(data['weight'])
        # control remark in the experiment
        remark = data['remark'].split()[0]
        controlValue = data['remark'].split()[1]
        controlUnit = data['remark'].split()[2]
        # raw data path and processed data path
        rawDataPath = os.path.join(E_MODUL_RAWDATA_PATH, nameOfExperiment)
        nameOfProcessedDataFile = 'processed_' + make_valid_uri_from_string(nameOfExperiment) + '.csv'
        processedDataPath = os.path.join(E_MODUL_RAWDATA_PATH, nameOfProcessedDataFile)

        g = My_world.as_rdflib_graph()

        with lebedigital_concrete:
            g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity('Experiment_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity.iri))
            )
            )
            # add specimens as instances in class Specimen
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('Specimen_' + make_valid_uri_from_string(specimenName)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen.iri))
                )
            )
            # add testers (prüfers) as instances in class Agent
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Agent('Agent_' + operator).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.Agent.iri))
                )
            )
            # add person/tester in class DesignativeName
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.DesignativeName('DesignativeName_' + operator).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.DesignativeName.iri))
                )
            )
            # add start date as instances in class Day
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Day('Day_' + make_valid_uri_from_string(str(operatorDay))).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.Day.iri))
                )
            )
            # add controls as instances in class ForceRate
            g.add(
                (
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate('ForceRate_' + make_valid_uri_from_string(nameOfExperiment) + remark).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate.iri))
                )
            )
            # add unit in MeasurementUnitOfForceRate
            g.add(
                (
                    URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + make_valid_uri_from_string(nameOfExperiment + controlUnit)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate.iri))
                )
            )
            # add length of the specimen in class Length
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Length('Length_' + make_valid_uri_from_string(specimenName)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.Length.iri))
                )
            )
            # add Diameter of the specimen in class Diameter
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Diameter('Diameter_' + make_valid_uri_from_string(specimenName)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.Diameter.iri))
                )
            )
            # add Weight of the specimen in class Mass
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Mass('Mass_' + make_valid_uri_from_string(specimenName)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.Mass.iri))
                )
            )
    #        
            # add specimen region as instances in class MeasurementRegion
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + make_valid_uri_from_string(specimenName)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.MeasurementRegion.iri))
                )
            )
            # add weight, diameter, length, dataset path, tester name value, force rate value, day
            # in class InformationBearingEntity
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenLength))).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenDiameter))).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenMass))).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + operator).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName) + controlValue ).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(str(operatorDay))).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )

            # add processed data into class bfo:BFO_0000015
            g.add(
                (
                    URIRef(urllib.parse.unquote(OBO.BFO_0000015('BFO_0000015_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(OBO.BFO_0000015.iri))
                )
            )
            # add processed data into class AnalysedDataSet
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet('AnalysedDataSet_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet.iri))
                )
            )
            # add processed data into class InformationBearingEntity
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet.iri))
                )
            )
            # adding object properties

            # the experiment has output raw data
            g.add(
                (
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity('Experiment_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_output.iri)), 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri))
                )
            )

            # experiment has operator
            g.add(
                (
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity('Experiment_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_agent.iri)), 
                    URIRef(urllib.parse.unquote(CCO.Agent('Agent_' + operator).iri))
                )
            )

            # specimen is the input of the experiment
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('Specimen_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.is_input_of.iri)), 
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity('Experiment_' + make_valid_uri_from_string(nameOfExperiment)).iri))
                )
            )

            # ForceRate is_input_of experiment
            g.add(
                (
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate('ForceRate_' + make_valid_uri_from_string(nameOfExperiment) + remark).iri)), 
                    URIRef(urllib.parse.unquote(CCO.is_input_of.iri)), 
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity('Experiment_' + make_valid_uri_from_string(nameOfExperiment)).iri))
                )
            )
            # Experiment occures_on Day
            g.add(
                (
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity('Experiment_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.occures_on.iri)), 
                    URIRef(urllib.parse.unquote(CCO.Day('Day_' + make_valid_uri_from_string(str(operatorDay))).iri))
                )
            )

            # specimen obo:BFO_0000051 MeasurementRegion
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('Specimen_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.BFO_0000051.iri)), 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + make_valid_uri_from_string(specimenName)).iri))
                )
            )
            # MeasurementRegion obo:BFO_0000086 Diameter, Length, Mass
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0000086.iri)), 
                    URIRef(urllib.parse.unquote(CCO.Diameter('Diameter_' + make_valid_uri_from_string(specimenName)).iri))
                )
            )

            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0000086.iri)), 
                    URIRef(urllib.parse.unquote(CCO.Diameter('Length_' + make_valid_uri_from_string(specimenName)).iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0000086.iri)), 
                    URIRef(urllib.parse.unquote(CCO.Diameter('Mass_' + make_valid_uri_from_string(specimenName)).iri))
                )
            )

            # Diameter, Length, Mass, ForceRate, RawDataSet, DesignativeName obo:RO_0010001 InformationBearingEntity
            # Day obo:RO_0010001 InformationBearingEntity
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Diameter('Diameter_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenDiameter))).iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Diameter('Length_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenLength))).iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Diameter('Mass_' + make_valid_uri_from_string(specimenName)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenMass))).iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate('ForceRate_' + make_valid_uri_from_string(nameOfExperiment) + remark).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + make_valid_uri_from_string(nameOfExperiment + controlUnit)).iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.DesignativeName('DesignativeName_' + operator).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)),
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + operator).iri))
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Day('Day_' + make_valid_uri_from_string(str(operatorDay))).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)),
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(str(operatorDay))).iri))
                )
            )
            # Agent designated_by DesignativeName
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.Agent('Agent_' + operator).iri)), 
                    URIRef(urllib.parse.unquote(CCO.designated_by.iri)),
                    URIRef(urllib.parse.unquote(CCO.DesignativeName('DesignativeName_' + operator).iri))
                )
            )
            # RawDataSet cco:is_input_of bfo:BFO_0000015
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri)), 
                    URIRef(urllib.parse.unquote(CCO.is_input_of.iri)), 
                    URIRef(urllib.parse.unquote(OBO.BFO_0000015('BFO_0000015_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri))
                )
            )
            # bfo:BFO_0000015 cco:has_output mseo:AnalysedDataSet
            g.add(
                (
                    URIRef(urllib.parse.unquote(OBO.BFO_0000015('BFO_0000015_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_output.iri)), 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet('AnalysedDataSet_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri))
                )
            )

            # mseo:AnalysedDataSet obo:RO_0010001 cco:InformationBearingEntity
            g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet('AnalysedDataSet_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri))
                )
            )

            # add data properties
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenMass))).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)),
                    Literal(specimenMass)
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenDiameter))).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)),
                    Literal(specimenDiameter)
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(specimenName + str(specimenLength))).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)),
                    Literal(specimenLength)
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(nameOfExperiment) + 'specimen_dat').iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_URI_value.iri)), 
                    Literal(rawDataPath)
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + operator).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_text_value.iri)), 
                    Literal(operator)
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + make_valid_uri_from_string(nameOfExperiment + controlUnit)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)), 
                    Literal(controlValue)
                )
            )
            g.add(
                (
                    URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + make_valid_uri_from_string(nameOfExperiment + controlUnit)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_text_value.iri)), 
                    Literal(controlUnit)
                )
            )
            # cco: InformationBearingEntity of processed dataset has filepath 
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + 'processed_' + make_valid_uri_from_string(nameOfExperiment)).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_URI_value.iri)), 
                    Literal(processedDataPath)
                )
            )
            # Day has_datetime_value
            g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + make_valid_uri_from_string(str(operatorDay))).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_datetime_value.iri)), 
                    Literal(operatorTime)
                )
            )

            g.serialize(destination=locationOfKnowledgeGraph, format="turtle")

    except:
        print('Please check the file path')

    


# metadata_ontology_mapping('C:\\Users\\vdo\\Desktop\\LeBeDigital\\Code\\minimum_working_example\\ModelCalibration\\usecases\\Concrete\\Example\\emodul\\metadata_yaml_files\\BA Los M V-4.yaml','C:\\Users\\vdo\\Desktop\\LeBeDigital\\Code\\minimum_working_example\\ModelCalibration\\usecases\\Concrete\\Example\\emodul\\triples\\emodul_knowledge_graph.ttl')