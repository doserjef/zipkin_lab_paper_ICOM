# Data

1. `guild-info.txt`: table that contains information on different structural, functional, and compositional bird guilds. This classification comes from [O'Connell et al. 2000](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/1051-0761%282000%29010%5B1706%3ABGAIOE%5D2.0.CO%3B2). 
2. `pa-data-bundle.rda`: eBird, BBS, and covariate data for running the ICOM to a community of interior forest obligate bird species. For initial model testing and fitting, this only contains data from PA. 
3. `spatial-covariates.rda`: covariate data that are saved separately from the `pa-data-bundle.rda` to avoid having to always re-calculate the variables. Just a temporary file, as all of the covariate data is eventually stored in `pa-data-bundle.rda`. 
