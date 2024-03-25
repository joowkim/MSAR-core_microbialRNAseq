# Metagenome Shotgun data QC/gene annotation pipeline

## Overview of the workflow

1. Rename fastq files. 
   - e.g. `GC.TN.1436_***.fastq.gz -> CRIS72_***.fastq.gz`
2. Run QC tools
   1. `fastqc` on each sample - raw fastq files 
   2. `fastp` to trim adapter sequences and low quality reads
      1. below options used for `fastp`
         - `--detect_adapter_for_pe`
         - `--qualified_quality_phred 20`
         - `--adapter_fasta $adapter`
   
3. `multiqc` for summarizing the outputs of qc tools
4. Removal of host genome using `bowtie2` 
   - i.e. align reads against host genome.
     1. below options used for `bowtie2` and `samtools`
     2. `bowtie2 -x [index_file] -1 R1 -2 R2 | samtools sort -o bam_file`
     3. `samtools view -b -f 12 -F 256 bam_files > both_unmapped.bam`
     4. `samtools fastq both_unmapped.bam -1 sample_host_remove.R1.fastq.gz -2 sample_host_remove.R2.fastq.gz`
     5. `sample_host_remove.R[12].fastq.gz` files are host removed fastq files
     6. [referece](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences) 
5. Obtain host genome filtered reads (`host_remove.R[12].fastq.gz`) / unfiltered reads (fastp output)
6. Genome assebmly using megahit using trimmed reads from raw fastq files - fastp output files. - Not `host_remove.fastq.gz` files
   - `final_contigs.fa` obtained.
   - `--k-min 27 --k-max 47` options used.
7. Gene annotation using prokka
    - Execute prokka - `final_contigs.fa` as input
    - default option used

## Visual of the workflow


## Dependency

`Nextflow`

`bowtie2`

`samtools`

`fastqc`

`fastp`

`multiqc`

`megahit`

`prokka`

`singularity`