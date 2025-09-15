#!/bin/bash

# Ethan Ross
# 15.01.25
# In this script I will download the taxonomy .sql file from NCBI into a scrath folder on Maxwell

# Here used decona on some test data using SLURM

#SBATCH --cpus-per-task=10 
#SBATCH --mem-per-cpu=5000
#SBATCH --partition=uoa-compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r01er21@abdn.ac.uk 

module load r/4.4.0

R

setwd("/uoa/scratch/users/r01er21/taxonomizr")

library(taxonomizr)

prepareDatabase("accessionTaxa.sql")
