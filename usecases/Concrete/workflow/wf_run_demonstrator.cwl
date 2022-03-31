class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
  
inputs:
  - id: mix_wz
    type: float
    doc: water cement ratio of the mix design
    default: 1
  - id: radius
    type: float
    doc: radius of test cylinder in mm
    default: 75
  - id: height
    type: float
    doc: height of test cylinder in mm
    default: 300
    
outputs:
  - id: KPI_CO2_emission
    outputSource:
      - structural_simulation/KPI_CO2_emission
    type: float
    doc: C02 emission of the global structure
  - id: KPI_demoulding_time
    outputSource:
      - structural_simulation/KPI_demoulding_time
    type: float
    doc: min time for demoulding
    
steps:
  - id: mix_design_performance_prediction_test
    in:
      - id: mix_wz
        source: mix_wz
      - id: radius
        source: radius
    out:
      - id: CO2
    run: mix_design_performance_prediction.cwl
  - id: structural_simulation
    in:
      - id: CO2
        source:
          - mix_design_performance_prediction/CO2          
      - id: height
        source: height
    out:
      - id: KPI_CO2_emission
      - id: KPI_demoulding_time
    run: structural_simulation.cwl
requirements: []
