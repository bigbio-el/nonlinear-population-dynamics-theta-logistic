# Nonlinear population dynamics of the Chukchi Sea polar bear population

MATLAB implementation of a reduced theta-logistic model for the Chukchi Sea polar bear (*Ursus maritimus*) subpopulation under habitat change and fixed harvesting.

## Overview

This repository accompanies a study of nonlinear population dynamics in the Chukchi Sea polar bear population. The model represents density-dependent population growth, a time-dependent decline in environmental carrying capacity, and fixed annual harvesting.

It reproduces the parameter calculations, deterministic scenario analyses, harvest-threshold calculations, sensitivity analyses, and figures reported in the associated manuscript.

The model is intended for transparent scenario analysis. It is not a population viability analysis, a calendar-dated forecast, or a basis for setting sustainable-harvest quotas.

## Contents

* `polar_bear_population_dynamics_theta.m` — main MATLAB script for reproducing the analyses and figures.

## Requirements

* MATLAB
* No external toolboxes are required.

## Running the code

1. Clone or download this repository.
2. Open MATLAB and set this repository as the current folder.
3. Run:

```matlab
polar_bear_population_dynamics_theta
```

Output figures are automatically written to the `figures/` directory.

## Data sources

No new empirical data were generated or collected. The model inputs are derived from published sources, principally:

* Regehr, E. V. et al. (2018). *Integrated population modeling provides the first empirical estimates of vital rates and abundance for polar bears in the Chukchi Sea*. **Scientific Reports**, 8, 16780. https://doi.org/10.1038/s41598-018-34824-7

* Regehr, E. V. et al. (2021). *Demographic risk assessment for a harvested species threatened by climate change: polar bears in the Chukchi Sea*. **Ecological Applications**, 31, e02461. https://doi.org/10.1002/eap.2461

All numerical inputs are stated explicitly in the MATLAB script and are traceable to the associated manuscript and its references.

## Authors

Evgeniia Lavrenteva and Michael Binns

## Citation

If you use this code, please cite the associated manuscript. 
