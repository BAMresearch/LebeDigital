class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - python3
inputs:
  - id: height
    type: float
    doc: >-
      CO2 emission per m³ of concrete (make sure to upscale from mortar to
      concrete)
  - id: radius
    type: float
    doc: >-
      CO2 emission per m³ of concrete (make sure to upscale from mortar to
      concrete)
  - default:
      class: File
      location: cylinder_simulation.py
    id: script
    type: File
    
outputs:
  - id: force
    doc: bla
    type: float
    outputBinding:
      loadContents: true
      glob: output_KPI_CO2.txt
      outputEval: |
        ${
            var output = self[0].contents;
            var KPI_CO2_string = output.split("KPI_CO2:")[1];
            KPI_CO2_string = KPI_CO2_string.split("\n")[0];
            return parseFloat(KPI_CO2_string);
        }
doc: |
  Macroscopic simulation of the precast element
arguments:
  - $(inputs.script)
  - '--CO2'
  - $(inputs.CO2)
  - '--output'
  - output_KPI_CO2.txt
requirements:
  - class: InlineJavascriptRequirement
