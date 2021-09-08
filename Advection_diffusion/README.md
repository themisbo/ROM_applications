# Advection - diffusion

## Data

We derive the data for the advection - diffusion example from [[1]](https://www.pnas.org/content/116/31/15344).
For the data creation we use the datadrivenpdes package, which requires an older version of tensorflow. Therefore, we need to create a new virtual environment, install the new requirements and the datadrivenpdes package:

```sh
conda create -n advection python=3.7
conda activate advection
pip install -r requirements.txt
pip install git+git//github.com/google-research/data-driven-pdes
```

Then, within the [Advection_diffusion folder](.) we can create the data by running:
```sh
python Advection_data_creation.py
```

## Pipeline

For the pipeline make sure that you have reverted back to the ROM virtual environment:
```sh
conda deactivate
conda activate ROM
```

The code can be run with the [jupyter notebook](advection_pipeline.ipynb).<br />
The metrics are saved into the [txt file](metrics.txt).<br />
The plots are saved into the [plots folder](plots/).


## References

[1] Bar-Sinai, Y., Hoyer, S., Hickey, J., & Brenner, M. P. (2019). Learning data-driven discretizations for partial differential equations. _Proceedings of the National Academy of Sciences_, _116_(31), 15344-15349.
