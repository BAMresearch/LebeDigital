import pytest
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

from lebedigital.demonstrator_calibration.forward_solvers import HomogenizationSolverWrapper, HydrationSolverWrapper

def test_hydration_solver_wrapper():
    # -- observed inputs
    inp_solver = {}
    inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = [0,5000,10000,20000,100000]

    # -- latents -----
    # b = np.array([2.916,2.4229,5.554,5])
    # std = np.array([1.9956, 247.6045,   1.8181,   2.5245]) 
    # mean = np.array([  2.8128, 124.1033,   3.4967,   3.6444])
    # b = (b-mean)/std
    
    b = np.array([  2.91, np.log(2.422e-03),   3.4967,   3.6444, 4.7002])


    hydration_solver = HydrationSolverWrapper()
    heat_list = hydration_solver.solve(latents=b,inp_solver=inp_solver)
    #heat_list = hydration_solver_wrapper(b,inp_solver)
    print(f'heat_list = {heat_list}')

    # -- expected outputs
    heat_list_exp =[0.0,         0.06231378, 0.1284133,  0.27287546, 2.28620259]
    # assert the values are approximately equal
    # write assert statement also
    assert np.allclose(heat_list,heat_list_exp,atol=1e-3), "The heat list is not equal to the expected values"

def test_homogenization_solver():
    latents = [30,3]
    homogenization_solver = HomogenizationSolverWrapper()
    result = homogenization_solver.solve(latents=latents)
    print(f'result = {result}')
    result_correct = [51082128028.566986, 38101522.84263957]
    assert np.allclose(result,result_correct,atol=1e-3), "The homogenization solver is not working properly"


#test_homogenization_solver()
#test_hydration_solver_wrapper()