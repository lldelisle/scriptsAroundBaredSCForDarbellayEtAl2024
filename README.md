# scriptsAroundBaredSCForDarbellayEtAl2024

All scripts used around baredSC (from baredSC run to figure) in Darbellay et al. 2024.

Tables with raw counts of Col2a1 and GFP can be found in the [input directory](./inputs/).

The [pipeline](./pipeline.sh) will create a conda environment, run baredSC on a slurm cluster. Then it will plot the results with the R scripts available in [this directory](./r_scripts). The result plots are in the [plots directory](./plots/).
