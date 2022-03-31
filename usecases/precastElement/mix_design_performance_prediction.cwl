#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  Prediction of mix design properties
baseCommand: ["python3"]
requirements:
  - class: InlineJavascriptRequirement

arguments: ["$(inputs.script)",
                "-wz", "$(inputs.mix_wz)",
                "-o", "outputCO2.txt"
            ]

inputs:
  script:
    type: File
    default:
      class: File
      location: mix_design_performance_prediction.py
  mix_wz:
    type: float
    doc: "water cement ratio of the mix"

outputs:
  CO2:
    type: float
    doc: "CO2 emission per m³ of concrete (make sure to upscale from mortar to concrete)"
    outputBinding:
      glob: "outputCO2.txt"
      loadContents: true
      outputEval: |
        ${
            var output = self[0].contents;
            var CO2_string = output.split("CO2:")[1];
            CO2_string = CO2_string.split("\n")[0];
            return parseFloat(CO2_string);
        }