This folder contains the R (and python) codes to replicate the results of our paper titled "Sequential Transport for Causal Mediation Analysis".

The scripts are in the `scripts/` folder.

- `01-1-gaussian-2mediators.R`: Toy example with 2 mediators (example with a single run of simulations).
- `01-2-gaussian-2mediators-simul.R`: Toy example with 2 mediators: Monte Carlo Simulations
- `02-gaussian-3mediators.R`: Toy example with 3 mediators: Monte Carlo simulations
- `03-compas.R`: Application to the COMPAS recidivism dataset

The results of estimations are saved in the `output/` folder, once the codes are run.
If you turn the boolean variables saving graphs to `TRUE`, the graphs are exported in the `figs/` folder.

# Reproduction of Figures and Tables

How to reproduce specific figures from the paper:

- Figure 2: `scripts/02-gaussian-3mediators.R`
- Figure 4: `scripts/03-compas.R.R`
- Figure 5 (Appendix): `scripts/01-1-gaussian-2mediators.R`
- Figure 6 (Appendix): `scripts/01-1-gaussian-2mediators.R`
- Figure 7 (Appendix): `scripts/01-1-gaussian-2mediators.R`
- Figure 8 (Appendix): `scripts/01-2-gaussian-2mediators-simul.R`
- Figure 9 (Appendix): `scripts/01-2-gaussian-2mediators-simul.R`
- Figure 10 (Appendix):  `scripts/01-2-gaussian-2mediators-simul.R`


How to reproduce specific tables from the paper:

- Table 1: `scripts/03-compas.R.R`
- Table 2 (Appendix): `scripts/01-1-gaussian-2mediators.R`