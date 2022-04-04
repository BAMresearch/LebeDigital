class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: mix_wz
outputs:
  - id: force
    outputSource:
      - cylinder_simulation/force
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
    
  - id: cylinder_simulation
    in:
      - id: height
        source: forward_model/height
      - id: radius
        source: forward_model/radius
    out:
      - id: force
    run: cylinder_simulation.cwl
requirements: []
