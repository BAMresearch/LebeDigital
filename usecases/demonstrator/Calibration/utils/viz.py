# Viz files needed:
# 1. Trace plot of for the samples from E step.
# 2. Evolution of paramters phi and the grads
# 3. probabilistc map between x and b plot.
# 4. P{rediction of b and posterior predictive of the solver output subsequently
import json
import os, sys
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc

plt.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath,bm}']
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")


class Viz:
    @staticmethod
    def posterior_predictive():
        pass


def plot_constraints_and_objective(objective_func: callable, constraints: list[callable], x_bounds: tuple[float, float],
                                   y_bounds: tuple[float, float], best_x: float = None, best_y: float = None,
                                   x_steps: int = 100,
                                   y_steps: int = 100,**kwargs) -> None:
    """
    Plots the objective function and constraints on a 2D plot. Limited to 2 design variables

    Parameters:
    -----------
    objective_func : callable
        The objective function to be optimized.
    constraints : list[callable]
        A list of constraint functions. Each constraint function should return a value greater than or equal to 0
        for valid input points (x, y).
    x_bounds : tuple[float, float]
        A tuple of the form (x_min, x_max) representing the bounds of the x-axis.
    y_bounds : tuple[float, float]
        A tuple of the form (y_min, y_max) representing the bounds of the y-axis.
    best_x : float
        The x-coordinate of the optimal point.
    best_y : float
        The y-coordinate of the optimal point.
    x_steps : int, optional
        The number of steps to take along the x-axis, default is 100.
    y_steps : int, optional
        The number of steps to take along the y-axis, default is 100.

    Returns:
    --------
    None
    """

    x_grid = np.linspace(x_bounds[0], x_bounds[1], x_steps)
    y_grid = np.linspace(y_bounds[0], y_bounds[1], y_steps)
    obj_vals = np.zeros((len(x_grid), len(y_grid)))
    con_vals = [np.zeros((len(x_grid), len(y_grid))) for _ in constraints]
    for i, x in enumerate(x_grid):
        for j, y in enumerate(y_grid):
            obj_vals[i, j] = objective_func(x, y,**kwargs)
            for k, constraint in enumerate(constraints):
                con_vals[k][i, j] = constraint(x, y,**kwargs)
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    obj_contour = ax.contourf(x_grid, y_grid, obj_vals, levels=20, cmap='inferno')
    for k, con_vals_k in enumerate(con_vals):
        ax.contour(
            x_grid,
            y_grid,
            con_vals_k,
            levels=[0], # this ensures that just a sharp deviding line is present, kind of like indicator function. marks where there is a change of sign
            colors=[f"C{k}"],
            linestyles=["dashed"],
            label=f"Constraint {k + 1}",
        )
    if best_x is not None:
        ax.scatter(best_x, best_y, s=100, marker='*', color='k')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Objective Function and Constraints')
    #ax.legend()
    plt.colorbar(obj_contour)
    plt.show()

    return fig, obj_vals, con_vals
