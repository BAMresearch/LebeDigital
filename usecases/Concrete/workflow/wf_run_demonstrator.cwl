class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: mix
    'sbg:x': 13.918392181396484
    'sbg:y': 156.72906494140625
  - id: probeye_module
    'sbg:x': 13.9183931350708
    'sbg:y': -151.86886596679688
  - id: simulation_module
    'sbg:x': 345.64007568359375
    'sbg:y': 374.63671875
  - id: forward_model_definition
    'sbg:x': 10.438794136047363
    'sbg:y': -3.2588558197021484
outputs:
  - id: optimal_model_parameters
    outputSource:
      - calibration/optimal_model_parameters
    type: null
    'sbg:x': 979.308837890625
    'sbg:y': -69.64078521728516
steps:
  - id: data_querry
    in:
      - id: mix
        source: mix
    out:
      - id: experimental_data
      - id: meta_data
    run: data_querry.cwl
    'sbg:x': 251.5
    'sbg:y': 153.171875
  - id: calibration
    in:
      - id: displ_list
        source: data_querry/experimental_data
      - id: force_list
        source: data_querry/experimental_data
      - id: forward_model_definition
        source: forward_model_definition
      - id: meta_data
        source: data_querry/meta_data
      - id: probeye_module
        source: probeye_module
      - id: sim_force_list
        source: forward_model/sim_force_list
    out:
      - id: displ_list
      - id: model_parameters
      - id: optimal_model_parameters
    run: calibration.cwl
    'sbg:x': 540.0018310546875
    'sbg:y': -32.56990051269531
  - id: forward_model
    in:
      - id: displ_list
        source: calibration/displ_list
      - id: force
        source: cylinder_simulation/force
      - id: model_parameters
        source: calibration/model_parameters
      - id: simulation_module
        source: simulation_module
    out:
      - id: displacement
      - id: model_parameters
      - id: sim_force_list
    run: forward_model.cwl
    'sbg:x': 630.1891479492188
    'sbg:y': 172.96719360351562
  - id: cylinder_simulation
    in:
      - id: height
        source: forward_model/model_parameters
      - id: radius
        source: forward_model/model_parameters
      - id: displacement
        source: forward_model/displacement
      - id: youngs_modulus
        source: forward_model/model_parameters
      - id: poissions_ratio
        source: forward_model/model_parameters
    out:
      - id: force
    run: cylinder_simulation.cwl
    'sbg:x': 789.8375244140625
    'sbg:y': 342.1076354980469
requirements: []
