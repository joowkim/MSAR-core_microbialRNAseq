#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=bwa_mem
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=168:00:00
#SBATCH -o bwa_mem.%A.o
#SBATCH -e bwa_mem.%A.e

module load nextflow/22.04.3
module load singularity/3.8.0

nextflow run bwa_mem.nf -c ./bwa_mem_conf/run.config -resume -profile slurm

module unload nextflow/22.04.3
module unload singularity/3.8.0