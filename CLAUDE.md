# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running Workflows

All workflows are run with Nextflow on an HPC cluster (SLURM). The general pattern is:

```bash
# Run custom gene mapping workflow
nextflow run custom_gene_mapping.nf -c custom_gene_mapping_conf/run.config -profile slurm

# Run QC workflow
nextflow run qc.nf -c qc_conf/run.config -profile slurm

# Run MetaPhlAn4 / metagenomics workflow
nextflow run metaphlan4.nf -c metaphlan4_conf/run.config -profile slurm
```

Resume a previous run with `-resume`. Override params at the command line with `--param_name value` (e.g., `--bwa_index /path/to/index`, `--run_kraken true`).

## Architecture

This repo contains multiple independent Nextflow DSL2 workflows for shotgun metagenomics. Each workflow has its own `.nf` file and a paired `*_conf/` directory.

### Workflows

| Workflow | File | Purpose |
|---|---|---|
| Custom gene mapping | `custom_gene_mapping.nf` | QC → trim → screen → BWA-MEM alignment → optional Kraken2 + kraken-biom |
| QC | `qc.nf` | FastQC + optional fastp + FastQ Screen (standalone, not modularized) |
| MetaPhlAn4 / HUMAnN | `metaphlan4.nf` | QC → trim → host removal (bowtie2) → MetaPhlAn4, HUMAnN3, Kraken2, MEGAHIT (feature-flagged) |
| BAM to Kraken | `from-bam-to-kraken.nf` | Starts from existing BAMs |

### Config Structure (per workflow)

Each `*_conf/` directory contains:
- `run.config` — default params (samplesheet path, feature flags, genome choice) + includes the others
- `cluster.config` — SLURM executor settings, Singularity, queue size/rate limits
- `processes.config` — resource labels (`process_low/medium/high/dual`, `memory_medium/high`) and container assignments
- `reference.config` — paths to reference indexes (bowtie2, bwa, etc.)

### Modules (`modules/`)

`custom_gene_mapping.nf` uses modularized processes (one process per file). `qc.nf` and `metaphlan4.nf` define processes inline.

### Samplesheet Format

All workflows read a CSV samplesheet (`samplesheet.csv` in `launchDir` by default):

```
sample,fq1,fq2
SampleA,/path/to/A_R1.fastq.gz,/path/to/A_R2.fastq.gz
SampleB,/path/to/B_R1.fastq.gz,   # empty fq2 = single-end
```

SE/PE is auto-detected: empty `fq2` → `single_end: true`. The `meta` map (`[sample_name:, single_end:]`) is passed through most channels.

### Output

All `publishDir` paths use `${launchDir}/analysis/<tool>/`. Outputs land next to where you launched Nextflow.

### Software

Most tools are loaded via HPC `module` directives. MultiQC uses a Singularity container (`processes.config`). Singularity is enabled globally in `cluster.config`.

`kraken_biom` is an exception — it uses a Python venv at `~/beegfs/python_env/kraken2/` invoked directly via its full binary path (no `module`, no `source activate`).

### Key Parameters (custom_gene_mapping)

| Param | Default | Description |
|---|---|---|
| `params.samplesheet` | `${launchDir}/samplesheet.csv` | Input CSV |
| `params.bwa_index` | `""` | Path to BWA index (required) |
| `params.run_name` | `"multiqc_report"` | MultiQC output filename |
| `params.host_genome` | `"mouse"` | Host genome key for fastq_screen |
| `params.run_kraken` | `false` | Enable filt_bam → get_mapped_reads → kraken2 → kraken_biom branch |
