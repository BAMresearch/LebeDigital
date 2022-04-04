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
      location: forward_model.py
  mix_wz:
    type: float
    doc: "water cement ratio of the mix"

outputs:
  height:
    type: float
  radius:
    type: float



