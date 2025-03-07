#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=nf-bwa_mem
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=168:00:00
#SBATCH -o slurm.%A.o
#SBATCH -e slurm.%A.e

module unload nextflow/23.10.0
module unload singularity/3.8.0
module unload java/21.0.6

nextflow run bwa_mem.nf -c ./bwa_mem_conf/run.config -resume -profile slurm

module unload nextflow/23.10.0
module unload singularity/3.8.0
module unload java/21.0.6