# Modelling the impact of poor-quality antimalarials on *Plasmodium falciparum* resistance dynamics

R code and selected model outputs associated with the manuscript **“Modelling the impact of poor-quality antimalarials on Plasmodium falciparum resistance dynamics”** by Aronrag Meeyai, Sean Cavany, Paul N. Newton and Lisa White.

## Overview

This repository contains a deterministic compartmental transmission model with interacting drug-sensitive and drug-resistant *P. falciparum* populations. The analyses explore how antimalarial active pharmaceutical ingredient (API) content, delays to treatment, transmission intensity and drug resistance influence malaria incidence and the emergence of resistance.

## Repository contents

- `nice_MQmalaria_Ci_Prec_gammaa_gammap_Function_A_B_C_D_E_Influence_withcomments_19Mar2019_forPlotHeatMap_recheck8April2019.R` — core transmission model and example simulations.
- `nice_MQmalaria_Ci_Prec_gammaa_gammap_Function_A_B_C_D_E_Influence_withcomments_19Mar2019_forPlotHeatMap_recheck8April2019_sweep.R` — parameter-sweep simulations.
- `sensitivity_analysis.R` and `sensitivity_analysis_qq.R` — global sensitivity analyses, including Latin hypercube sampling and an API-content sweep.
- `process_heatmap_data.R`, `plot_heatmap.R`, `plot_heatmap2.R`, and other plotting scripts — preparation and visualisation of simulation outputs.
- `HeatMap_data.csv` and `tt_*.csv` files — selected model-output data.

## Software and use

The analyses were written in **R**. Dependencies used across the scripts include `deSolve`, `lhs`, `sensitivity`, `pbapply`, `ppcor`, `tidyverse`, `viridis`, `RColorBrewer`, `patchwork`, and `cowplot`; some scripts also use `parallel` and `grid`.

Scripts can be opened and run in R or RStudio after installing their required packages. Some scripts assume particular working directories, input files or previously created objects; check file paths and script-specific settings before running. The repository is an archive of research scripts, rather than a packaged one-command workflow.

## Publication and citation

Meeyai A, Cavany S, Newton PN, White L. *Modelling the impact of poor-quality antimalarials on Plasmodium falciparum resistance dynamics*.

For a fixed, citable version of this code, use the Zenodo DOI associated with the relevant GitHub release: **[add Zenodo DOI after archiving]**.

## Licence

This code is made available under the MIT License.
See the [LICENSE](LICENSE) file for details.