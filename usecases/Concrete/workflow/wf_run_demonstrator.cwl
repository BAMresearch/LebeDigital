class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: mix_wz
    type: float
    doc: water cement ratio of the mix design
    default: 1
outputs:
  - id: KPI_CO2_emission
    outputSource:
      - structural_simulation/KPI_CO2_emission
    type: float
    doc: C02 emission of the global structure
steps:
  - id: mix_design_performance_prediction
    in:
      - id: mix_wz
        source: mix_wz
    out:
      - id: CO2
    run: mix_design_performance_prediction.cwl
  - id: structural_simulation
    in:
      - id: CO2
        source:
          - mix_design_performance_prediction/CO2
    out:
      - id: KPI_CO2_emission
    run: structural_simulation.cwl
requirements: []
