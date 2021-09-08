import numpy as np
import matplotlib.pyplot as plt
from simulator_2d import run2d
import xarray as xr

num_of_samples = 20
dif_vector = np.linspace(0.05, 0.15, num_of_samples)

for i in range(num_of_samples):
    print(i)
    concentration_simulation, _, _, _ = run2d(dif=dif_vector[i])
    np.save(
        "data/zvals_dif_" + str(round(dif_vector[i], 3)) + ".npy",
        concentration_simulation,
    )

# Construct data array object and plot indicative frames
concentration_simulation, x, y, times = run2d(dif=dif_vector[i], vel=0.5)

conc_true = xr.DataArray(
    concentration_simulation,
    dims=["time", "x", "y"],
    coords={"time": times, "x": x, "y": y},
)

frames_check = [0, 25, 49]
conc_true[frames_check, :, :].plot(col="time", robust=True, aspect=1)
plt.show()
