#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --output=nextflow.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=8G

# activate conda environment for nextflow
eval "$(conda shell.bash hook)"
conda activate /home/dzhou/miniconda3/envs/nextflow

# cd into directory with .nf file
cd ~/scRNA_analysis

# export NXF_SINGULARITY_CACHEDIR='/home/dzhou/scRNA_analysis/envs/'

# run nextflow pipeline script
nextflow run main.nf -resume
