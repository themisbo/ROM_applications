# Polymer precipitaiton

## Data

We derive the data for the polymer precipitation example from [[1]](https://arxiv.org/abs/2007.07276).
For the data creation we use the dolfin package, which requires a new virtual environment. We create the new environment, activate it and install the requirements:

```sh
conda create -n precipitation python=3.9
conda activate precipitation
pip install -r requirements.txt
```

Then, within the [Polymer precipitation folder](.) we can create the data by running:
```sh
python create_data1.py
python create_data2.py
```

## Pipeline

For the pipeline make sure that you have reverted back to the ROM virtual environment:
```sh
conda deactivate
conda activate ROM
```

The code can be run with the [jupyter notebook](precipitation_pipeline.ipynb).<br />
The metrics are saved into the [txt file](metrics.txt).<br />
The plots are saved into the [plots folder](plots/).

## References

[1] Inguva, P., Mason, L. R., Pan, I., Hengardi, M., and Matar, O. K. (2020). Numerical simulation, clustering, and prediction of multicomponent polymer precipitation. Data-Centric Engineering, 1.
