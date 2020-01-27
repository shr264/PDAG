#!/bin/bash
#SBATCH --job-name=partial_sim_
#SBATCH -o=partial_sim_stats_%A_%a.out
#SBATCH -e=partial_sim_stats_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shr264@ufl.edu
#SBATCH --account=informatics_inst
#SBATCH --qos=informatics_inst
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2000mb
#SBATCH -t 240:00:00
[[ -d $SLURM_SUBMIT_DIR ]] & cd $SLURM_SUBMIT_DIR
date;hostname;pwd

module purge
ml gcc/6.3.0 R/3.4.0

Rscript simulations.R

date
