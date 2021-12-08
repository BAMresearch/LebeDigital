from xml.dom.minidom import NamedNodeMap
import rdflib
from SPARQLWrapper import SPARQLWrapper
from pathlib import Path
import os
import sys
from compression_mapping import german_to_english

baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
triplePath = os.path.join(baseDir0,'compression-processed-data/compression_Graph.ttl')
demoPath = os.path.join(baseDir0,'compression-processed-data/demo.xml')

if sys.platform == 'win32':
    prefixPath = 'file:///' + os.path.join(baseDir0,'compression-processed-data').replace('\\','/') + '/'
else:
    prefixPath = 'file://' + os.path.join(baseDir0,'compression-processed-data') + '/'

graph = rdflib.Graph()
graph.parse(triplePath, format='n3')

def input_compression_data_for_calibration(nameOfExperiment):
    nameOfExperiment = 'compression_experiment_' + german_to_english(nameOfExperiment.replace(' ','_').replace('.','_'))
    q1 = f"""
            prefix bwmd: <https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#>
            prefix mseo: <https://purl.matolab.org/mseo/mid/>
            prefix cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
            prefix obo: <http://purl.obolibrary.org/obo/>
            prefix con: <https://github.com/BAMresearch/ModelCalibration/blob/Datasets/usecases/Concrete/ConcreteOntology/Concrete_Ontology_MSEO.owl#>

            select ?rawdatapath
            where {{
                {{
                    select ?info
                    where {{
                        {{
                            select ?analyseddata
                            where{{
                                {{
                                    select ?bfo
                                    where {{
                                        {{
                                            select ?rawdata
                                            where {{
                                                mseo:{nameOfExperiment}
                                                cco:has_output
                                                ?rawdata
                                            }}
                                        }}
                                        ?rawdata
                                        cco:is_input_of
                                        ?bfo
                                    }}
                                }}
                                ?bfo
                                cco:has_output
                                ?analyseddata
                            }}
                        }}
                        ?analyseddata
                        obo:RO_0010001
                        ?info
                    }}
                    
                }}
                ?info
                cco:has_URI_value
                ?rawdatapath
            }}
            limit 1
        """
    results = graph.query(q1)
    processedDataPath = ''
    for result in results:
        if sys.platform == 'win32':
            processedDataPath = result['rawdatapath'].value
        else:
            processedDataPath = str(result['rawdatapath'])


    
    specimenParameterNames = ['Mass', 'Diameter','Length']
    specimenParameters = []
    for parameter in specimenParameterNames:
        q2 = f"""
            prefix bwmd: <https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#>
            prefix mseo: <https://purl.matolab.org/mseo/mid/>
            prefix cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
            prefix obo: <http://purl.obolibrary.org/obo/>
            prefix con: <https://github.com/BAMresearch/ModelCalibration/blob/Datasets/usecases/Concrete/ConcreteOntology/Concrete_Ontology_MSEO.owl#>

            select ?parametervalue
            where {{
                {{
                    select ?info
                    where {{
                        {{
                            select ?parameterclass
                            where {{
                                {{
                                    select ?parameterclass
                                    where {{
                                        {{
                                            select ?measurementregion
                                            where {{
                                                {{
                                                    select ?specimen
                                                    where {{
                                                        ?specimen
                                                        cco:is_input_of
                                                        mseo:{nameOfExperiment}
                                                    }}
                                                }}
                                                ?specimen
                                                obo:BFO_0000051
                                                ?measurementregion
                                            }}
                                        }}
                                        ?measurementregion
                                        obo:RO_0000086
                                        ?parameterclass
                                    }}
                                }}
                                ?parameterclass
                                a
                                cco:{parameter}
                            }}
                        }}
                        ?parameterclass
                        obo:RO_0010001
                        ?info
                    }}
                }}
                ?info
                cco:has_decimal_value
                ?parametervalue
            }}
            """
        results = graph.query(q2)
        for result in results:
            if sys.platform == 'win32':
                specimenParameters.append(result['parametervalue'].value)
            else:
                specimenParameters.append(float(str(result['parametervalue'])))
    return {
        'processedDataPath': processedDataPath,
        'specimenMass': specimenParameters[0],
        'specimenDiameter': specimenParameters[1],
        'specimenHeight': specimenParameters[2]
    }
def turtle_file_for_demo(nameOfExperiment):
    nameOfExperiment = 'compression_experiment_' + german_to_english(nameOfExperiment.replace(' ','_').replace('.','_'))
    q = f"""
            prefix bwmd: <https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#>
            prefix mseo: <https://purl.matolab.org/mseo/mid/>
            prefix cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
            prefix obo: <http://purl.obolibrary.org/obo/>
            prefix con: <https://github.com/BAMresearch/ModelCalibration/blob/Datasets/usecases/Concrete/ConcreteOntology/Concrete_Ontology_MSEO.owl#>
            select ?s ?p ?o
            where {{
                {{
                    ?s
                    ?p
                    mseo:{nameOfExperiment}
                }}
                union
                {{
                    mseo:{nameOfExperiment}
                    ?p
                    ?o
                }}
            }}
        """
    results = graph.query(q)
    results.serialize(destination=demoPath, format="xml")
turtle_file_for_demo('Hüsken Probe 1-2')