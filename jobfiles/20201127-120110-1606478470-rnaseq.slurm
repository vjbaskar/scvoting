#!/bin/bash
#SBATCH -p clincloud-express
#SBATCH -A gottgens-ccld-sl2-cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --job-name rnaseq
#SBATCH --output .slurm/20201127-120110-1606478470-rnaseq.log
#SBATCH --error .slurm/20201127-120110-1606478470-rnaseq.log
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user vs401
#SBATCH --time 12:00:00
#SBATCH -p clincloud-express
bash batch_scripts/jl_scrnaseq0.8.sh