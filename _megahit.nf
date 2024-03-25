nextflow.enable.dsl=2

process fastp {
    debug true
    tag "${sample_name}"
    // label "universal"
    cpus 8
    memory '8 GB'

    publishDir "${launchDir}/analysis/fastp/"

    module 'fastp/0.21.0'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}_trimmed.R{1,2}.fq.gz"), emit: trim_reads
    path("${sample_name}.fastp.json"), emit: json

    script:
    adapter = "/mnt/beegfs/kimj32/reference/adapters.fa"
    """
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    --thread ${task.cpus} \
    --qualified_quality_phred 20 \
    -o ${sample_name}_trimmed.R1.fq.gz \
    -O ${sample_name}_trimmed.R2.fq.gz \
    --adapter_fasta $adapter \
    --json ${sample_name}.fastp.json
    """
}


process fastqc {
    debug true
    tag "Fastqc on ${sample_name}"
    //label "universal" // getting cpu and memory usage from salmon.config - called universal
    cpus 8
    memory '8 GB'

    publishDir "${launchDir}/analysis/fastqc/", mode: "copy"

    module 'FastQC/0.11.9'

    input:
    tuple val(sample_name), path(reads)

    output:
    path ("*.zip"), emit: zips
    path ("*.html"), emit: htmls

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}


process multiqc {
    debug true
    tag "Multiqc on this project"

    cpus 2
    memory '2 GB'

    publishDir "${launchDir}/analysis/multiqc/", mode : "copy"

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    config_yaml = "/home/kimj32/config_defaults.yaml"
    """
    multiqc ${files} --filename "multiqc_report.megahit.html" --config ${config_yaml}
    """
}


process fastq_screen {
    debug true
    tag "Fastq-screen on ${sample_name}"

    cpus 8
    memory '16 GB'

    publishDir "${launchDir}/analysis/fastq_screen", mode : "copy"

    module 'FastQScreen/0.14.1'
    module 'bowtie2/2.3.4.1'

    input:
    tuple val(sample_name), path(reads)

    output:
    path("*.html")
    path("*.txt"), emit: fastq_screen_out

    // threads option is already defined in fastq_screeN_conf
    script:
    conf = "/mnt/beegfs/kimj32/polymerase/polymeraseDependencies/FastQ_Screen_Genomes/fastq_screen.conf"
    """
    fastq_screen --aligner bowtie2 \
    --conf ${conf} \
    ${reads[0]} \
    --outdir ./  \
    --threads ${task.cpus}
    """
}


process bowtie2 {
    debug true
    tag "bowtie2 on ${sample_name}"

    cpus 8
    memory '16 GB'

    publishDir "${launchDir}/analysis/bowtie2"

    module 'bowtie2/2.3.4.1'
    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}.mapped_unmapped.bam"), emit: bowtie2_mapped_unmapped_bam
    tuple val(sample_name), path("${sample_name}.both_unmapped.bam"), emit: bowtie2_bam_both_unmapped_bam
    path("${sample_name}.mapped_unmapped.stats"), emit: samtools_stats

    script:
    index = "/mnt/beegfs/kimj32/reference/mouse/gencode/GRCm38.p6/indexes/bowtie2/mouse"
    """
    bowtie2 -p ${task.cpus} -x ${index} \
    -1 ${reads[0]} -2 ${reads[1]} \
    | samtools sort -@ 4 -O BAM -o ${sample_name}.mapped_unmapped.bam

    samtools index -@ 4 ${sample_name}.mapped_unmapped.bam
    samtools stats -@ 4 ${sample_name}.mapped_unmapped.bam > ${sample_name}.mapped_unmapped.stats

    samtools view -@ 4 -b -f 12 -F 256 \
    ${sample_name}.mapped_unmapped.bam > ${sample_name}.both_unmapped.bam

    """
}


process samtools_stats {
    debug true
    tag "samtools on ${both_unmapped_bam}"

    cpus 8
    memory '8 GB'

    publishDir "${launchDir}/analysis/samtools_stats", mode: "copy"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(both_unmapped_bam)

    output:
    path("${sample_name}.both_unmapped.stats"), emit: samtools_stats
    tuple val(sample_name), path("${sample_name}.both_unmapped.sorted.bam"), emit: both_unmapped_bam

    script:
    """
    samtools sort -@ 4 -m 5G ${both_unmapped_bam} -O BAM -o ${sample_name}.both_unmapped.sorted.bam
    samtools index -@ 4 ${sample_name}.both_unmapped.sorted.bam
    samtools stats -@ 4 ${sample_name}.both_unmapped.sorted.bam > ${sample_name}.both_unmapped.stats
    """
}


process split_reads_from_unmapped {
    debug true
    tag "split reads - ${sample_name}"

    cpus 8
    memory '8 GB'

    publishDir "${launchDir}/analysis/split_reads", mode: "copy"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam_file)

    output:
    tuple val(sample_name), path("${sample_name}.host_remove.R{1,2}.fastq.gz"), emit: split_reads

    script:
    """
    samtools sort -n -m 5G -@ 2 ${bam_file} -o ${sample_name}.sorted.bam
    samtools fastq -@ 8 ${sample_name}.sorted.bam \
        -1 ${sample_name}.host_remove.R1.fastq.gz \
        -2 ${sample_name}.host_remove.R2.fastq.gz \
        -0 /dev/null -s /dev/null -n
    """
}


process megahit {
    debug true
    tag "${sample_name}"

    cpus 12
    memory '32 GB'

    publishDir (
        path: "${launchDir}/analysis/megahit",
        pattern: "**"
        saveAS: { filename -> filename.endsWith("final.contigs.fa") ? "${sample_name}/${sample_name}.contigs.fa" : filename}
    )

    module "megahit/1.2.9"

    input:
    tuple val(sample_name), path(reads)

    output:
    path("${sample_name}", type: "dir")
    tuple val(sample_name), path("${sample_name}/final.contigs.fa"), emit: final_contigs_fa

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} \
        -o ${sample_name} \
        -t ${task.cpus} \
        --k-min 27 \
        --k-max 47
    """
}


process prokka {
    debug true
    tag "${sample_name}"

    cpus 8
    memory '16 GB'

    publishDir "${launchDir}/analysis/prokka", mode: "copy"

    // use singluarity to avoid install required additional tools for prokka - e.g. cpan...

    input:
    tuple val(sample_name), path(fasta)

    output:
    path("*")

    script:
    """
    prokka ${fasta} \
    --outdir ${sample_name} \
    --prefix ${sample_name} \
    --cpus ${task.cpus}
    """
}


workflow {
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    fastqc(ch_reads)

    if (params.run_fastp) {
        fastp(ch_reads)
        fastq_screen(fastp.out.trim_reads)
        bowtie2(fastp.out.trim_reads)
        samtools_stats(bowtie2.out.bowtie2_bam_both_unmapped_bam)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.json, bowtie2.out.samtools_stats, samtools_stats.out.samtools_stats).collect() )
        split_reads_from_unmapped( samtools_stats.out.both_unmapped_bam )
        megahit(fastp.out.trim_reads)
        prokka(megahit.out.final_contigs_fa)

    } else {
        fastq_screen(ch_reads)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips).collect() )
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}