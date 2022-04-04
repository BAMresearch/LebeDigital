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
  - id: KPI_CO2_emission
    outputSource:
      - structural_simulation/KPI_CO2_emission
    type: float
    doc: C02 emission of the global structure
    'sbg:x': 571.2452392578125
    'sbg:y': 7
steps:
  - id: mix_design_performance_prediction
    in:
      - id: mix_wz
        source: mix_wz
    out:
      - id: CO2
    run: mix_design_performance_prediction.cwl
    'sbg:x': 124.109375
    'sbg:y': 0
  - id: structural_simulation
    in:
      - id: CO2
        linkMerge: merge_flattened
        source:
          - mix_design_performance_prediction/CO2
    out:
      - id: KPI_CO2_emission
    run: structural_simulation.cwl
    'sbg:x': 335.4375
    'sbg:y': 0
requirements: []
