#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=QC
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=12:00:00
#SBATCH -o QC.%A.o
#SBATCH -e QC.%A.e

module load nextflow/22.04.3
module load singularity/3.8.0

nextflow run qc.nf -c ./qc_conf/run.config -resume

module unload nextflow/22.04.3
module unload singularity/3.8.0