#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  Macroscopic simulation of the precast element
baseCommand: ["python3"]

requirements:
  - class: InlineJavascriptRequirement

arguments: ["$(inputs.script)",
                "--CO2", "$(inputs.CO2)",
                "--output", "output_KPI_CO2.txt"
            ]
inputs:
  script:
    type: File
    default:
      class: File
      location: structural_simulation.py
  CO2:
    type: float
    doc: "CO2 emission per m³ of concrete (make sure to upscale from mortar to concrete)"

outputs:
  KPI_CO2_emission:
    type: float
    doc: "C02 emission of the global structure"
    outputBinding:
      glob: "output_KPI_CO2.txt"
      loadContents: true
      outputEval: |
        ${
            var output = self[0].contents;
            var KPI_CO2_string = output.split("KPI_CO2:")[1];
            KPI_CO2_string = KPI_CO2_string.split("\n")[0];
            return parseFloat(KPI_CO2_string);
        }