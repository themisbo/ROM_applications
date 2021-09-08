# Falling films

## Data

The data for the falling films application are derived from [[1]](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/wave-patterns-in-film-flows-modelling-and-threedimensional-waves/1134E94525E4968C89DA31F8C8605E2B), and were simulated using a julia language port of wavemeker [[2]](https://www.sciencedirect.com/science/article/pii/S235271101830075X).<br />
In order to create the data we install the julia dependencies and run the script:

```sh
julia requirements.jl
julia prod_data.jl
```

## Pipeline

For the pipeline make sure that you have reverted back to the ROM virtual environment:
```sh
conda deactivate
conda activate ROM
```

The code can be run with the [jupyter notebook](films_pipeline.ipynb).<br />
The metrics are saved into the [txt file](metrics.txt).<br />
The plots are saved into the [plots folder](plots/).


## References

[1] Scheid, B., Ruyer-Quil, C., and Manneville, P. (2006).  Wave patterns in film flows: modelling and three-dimensional waves. Journal of Fluid Mechanics, 562:183–222.<br />
[2] Rohlfs, W., Rietz, M., and Scheid, B. (2018). Wavemaker: The three-dimensional wave simulation tool forfalling liquid films.SoftwareX, 7:211–216.