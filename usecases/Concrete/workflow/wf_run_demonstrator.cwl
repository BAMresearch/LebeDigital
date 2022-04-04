class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: mix_wz
    type: float
    doc: water cement ratio of the mix design
    default: 1
    'sbg:x': -183
    'sbg:y': -153
outputs:
  - id: force
    outputSource:
      - cylinder_simulation/force
    type: float
    doc: C02 emission of the global structure
    'sbg:x': 571.2452392578125
    'sbg:y': 7
steps:
  - id: forward_model
    in:
      - id: mix_wz
        source: mix_wz
    out:
      - id: height
    out:
      - id: radius
    run: forward_model.cwl
    'sbg:x': 124.109375
    'sbg:y': 0
  - id: cylinder_simulation
    in:
      - id: height
        source: forward_model/height
      - id: radius
        source: forward_model/radius
    out:
      - id: force
    run: cylinder_simulation.cwl
    'sbg:x': 335.4375
    'sbg:y': 0
requirements: []
