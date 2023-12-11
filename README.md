# On the Landau gauge ghost-gluon-vertex close to and in the conformal window
This repository contains the code used to prepare the plots and results included in [On the Landau gauge ghost-gluon-vertex close to and in the conformal window]().

## Instructions: Rerunning the solver for the Dyson-Schwinger-equations
- Install required dependencies (see below)
- Download the data from the [Zenodo data release]() and place it in `input`
- Run the analysis using `julia main.jl`

## Instructions: Recreating the plots
- Install required dependencies (see below)
- Download the data from the [Zenodo data release]() and place it in `input`
- Generate the plots by running `julia plots.jl`

## Plots

The plots are made using [Plots.jl](https://zenodo.org/record/7994271) via the [PGFPlotsX](https://github.com/KristofferC/PGFPlotsX.jl) backend which requires a LaTeX installation with the PGFPlots package.

## Requirements
- julia 1.9
- LaTeX (including PGFPlots)

## Acknowledgements
FZ has been supported in part by the Austrian Science Fund research teams grant STRONG-DM (FG1) and the STFC Grant No. ST/X000648/1.