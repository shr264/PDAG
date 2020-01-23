#!/bin/bash
#SBATCH --job-name=rstudio
#SBATCH --output=rstudio_%j.log
#SBATCH --partition=gui

#SBATCH --account=presnell
#SBATCH --qos=presnell

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
date; hostname; pwd;

unset XDG_RUNTIME_DIR

module purge; module load gui/2 rstudio

/apps/gui/2.0.0/bin/start_gui_app rstudio

date