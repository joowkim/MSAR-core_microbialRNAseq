# Microbial RNAseq using `megahit, metaphlan4` and `humann3` pipeline

## Overview of the workflow

1. Create a samplesheet file to execute the pipeline which should be a csv file following the format below:

| sample  | fq1                 | fq2                 |
| ------- | ------------------- | ------------------- |
| sampleA | sampleA_R1.fastq.gz | sampleA_R2.fastq.gz |

2. Run QC tools
   1. `fastqc` on each sample - raw fastq files 
   2. `fastp` to trim adapter sequences and low quality reads
      1. below options used for `fastp`
         - `--qualified_quality_phred 20`
         - `--adapter_fasta $adapter`
   3. `fastq_screen` on R1 fastq files to detect possible contaminants.   
3. `multiqc` for summarizing the output files of the qc tools
4. Removal of host genome using `bowtie2` 
   - i.e. align reads against host genome.
     1. below options used for `bowtie2` and `samtools`
     2. `bowtie2 -x [index_file] -1 R1 -2 R2 | samtools sort -o bam_file`
     3. `samtools view -b -f 12 -F 256 bam_files > both_unmapped.bam`
     4. `samtools fastq both_unmapped.bam -1 sample_host_remove.R1.fastq.gz -2 sample_host_remove.R2.fastq.gz`
     5. `sample_host_remove.R[12].fastq.gz` files are host removed fastq files
     6. [referece](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences) 
5. Obtain host genome filtered reads (`host_remove.R[12].fastq.gz`) / unfiltered reads (fastp output)
6. Genome assembly using `megahit` using the host sequences filtered reads with `--k-min 27 --k-max 47` option.
7. microbe classification using `metaphlan4` on host genome filtered reads with `-t rel_ab_w_read_stats` option.
8. concatenate host sequence filtered reads fastq files for `humann3` 
   - i.e.`cat sample1_R1.fastq.gz sample1_R2.fastq.gz > sample1.concat.fastq.gz`
   - see [reference](https://forum.biobakery.org/t/humann3-paired-end-reads/862)
9. Run `humann3` on concatenated fastq files


## How to execute the pipeline
Adjust the configureation files such as `metaphlan4_conf/run.config and cluster.config` After that, 
```
sbatch run_metaphlan4.slurm
```

## Dependency

`slurm` `Nextflow` `bowtie2` `samtools` `fastqc` `fastp` `multiqc` `megahit` `metaphlan4` `singularity`