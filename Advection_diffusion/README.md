# Advection - diffusion

## Data
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
