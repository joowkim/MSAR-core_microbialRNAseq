#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=nf-metagenome
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=168:00:00
#SBATCH -o slurm.%A.o
#SBATCH -e slurm.%A.e

module load nextflow/23.0.4.1
module load singularity/3.8.0

nextflow run metaphlan4.nf -c ./metaphlan4_conf/run.config -resume -profile slurm

module unload nextflow/23.0.4.1
module unload singularity/3.8.0