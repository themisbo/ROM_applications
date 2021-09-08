#!/usr/bin/env burgers

from matplotlib import pyplot as plt
import numpy as np
from datadrivenpdes.core import grids
import datadrivenpdes as pde
import tensorflow as tf


tf.enable_eager_execution()


def run2d(SIZE=64, dif=0.05, vel=0.5, start_x=5, start_y=5):

    equation = pde.advection.equations.FiniteVolumeAdvectionDiffusion(
        diffusion_coefficient=dif
    )
    grid = grids.Grid.from_period(size=SIZE, length=2 * np.pi)
    x, y = grid.get_mesh()

    initial_state = equation.random_state(grid, seed=7109179)

    time_step = equation.get_time_step(grid)
    times = time_step * np.arange(0, 200, 4)

    results = pde.core.integrate.integrate_times(
        model=pde.core.models.FiniteDifferenceModel(equation, grid),
        state=initial_state,
        times=times,
        axis=0,
    )

    conc = results["concentration"].numpy()

    return conc, x[:, 0], y[0], times


if __name__ == "__main__":
    conc[0, :, :].plot()
    conc[::49].plot(col="time", robust=True, aspect=1)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.gca()
    ax.pcolor(y, x, conc[-20, :, :])
    ax.set_title("final concentration")
    ax.axvline(4.0, color="red", linewidth=2)
    ax.axhline(4.0, color="red", linewidth=2)
