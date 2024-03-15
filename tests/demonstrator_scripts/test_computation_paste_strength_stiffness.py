import pytest

from lebedigital.demonstrator_scripts.computation_paste_strength_stiffness import computation_paste_strength_stiffness
from lebedigital.unit_registry import ureg

def test_computation_hydration_parameters():
    # load the files
    NN_path = 'demonstrator_scripts/input_for_tests/NN_model_homogenization_final.pt'
    cov_path = 'demonstrator_scripts/input_for_tests/cov_parameters_homogenization_final.csv'

    
    E,fc = computation_paste_strength_stiffness(slag_ratio=0.2,gaussian_mean=NN_path,
                                                                        gaussian_cov=cov_path,seed=5)
    
    print(f'The paste parameters are: E = {E}, fc = {fc}')
    assert E.magnitude == pytest.approx(10.3219, abs=1e-3)
    assert fc.magnitude == pytest.approx(32.5727, abs=1e-3)

if __name__ == "__main__":
    test_computation_hydration_parameters()
