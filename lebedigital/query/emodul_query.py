import rdflib
from pathlib import Path
import os
import sys

BASEDIR1 = Path(__file__).resolve().parents[1]
BASEDIR2 = Path(__file__).resolve().parents[2]
E_MODUL_KNOWLEDGE_GRAPH_PATH = os.path.join(os.path.join(os.path.join(os.path.join(BASEDIR2,'Example'),'emodul'),'triples'),'emodul_knowledge_graph.ttl')

graph = rdflib.Graph()
graph.parse(E_MODUL_KNOWLEDGE_GRAPH_PATH, format='turtle')

def input_emodul_data_for_calibration(nameOfExperiment):
    nameOfExperiment = 'Experiment_' + nameOfExperiment.replace(' ','_').replace('.','_')
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
        processedDataPath = result['rawdatapath'].value

    
    specimenParameterNames = ['Mass', 'Diameter', 'Length']
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
        'specimenLength': specimenParameters[2]
    }


print(input_emodul_data_for_calibration('BA Los M V-4'))