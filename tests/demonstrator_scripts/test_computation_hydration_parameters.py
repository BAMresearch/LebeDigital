import pytest
import os
from pathlib import Path
from lebedigital.demonstrator_scripts.computation_hydration_parameters import computation_hydration_parameters
from lebedigital.unit_registry import ureg

def test_computation_hydration_parameters():
    # load the files
    cwd = Path(os.getcwd())
    #NN_path = cwd / Path("tests/demonstrator_scripts/input_for_tests/NN_model_hydration_final.pt")
    #cov_path = cwd / Path("tests/demonstrator_scripts/input_for_tests/cov_parameters_hydration_final.csv")
    NN_path = 'demonstrator_scripts/input_for_tests/NN_model_hydration_final.pt'
    cov_path = 'demonstrator_scripts/input_for_tests/cov_parameters_hydration_final.csv'


    B1, B2, eta, E_act, Q_pot, T_ref = computation_hydration_parameters(slag_ratio=0.2,gaussian_mean=NN_path,
                                                                        gaussian_cov=cov_path,seed=5)

    print(f'The hydration parameters are: B1 = {B1}, B2 = {B2}, eta = {eta}, E_act = {E_act}, Q_pot = {Q_pot}, T_ref = {T_ref}')
    assert B1.magnitude == pytest.approx(0.000453555218)
    assert B2.magnitude == pytest.approx(2.97396597e-06)
    assert eta.magnitude == pytest.approx(4.245238)
    assert E_act.magnitude == pytest.approx(5653 * 8.3145)
    assert Q_pot.magnitude == pytest.approx(469243.698)
    assert T_ref.magnitude == pytest.approx(20)

