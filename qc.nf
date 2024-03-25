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
    --trim_poly_g \
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

    cpus 1
    memory '1 GB'

    publishDir "${launchDir}/analysis/multiqc/", mode : "copy"

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    config_yaml = "/home/kimj32/config_defaults.yaml"
    """
    multiqc ${files} --filename "multiqc_report.qc.html" --config ${config_yaml}
    """
}


process fastq_screen {
    debug true
    tag "Fastq-screen on ${sample_name}"

    cpus 4
    memory '8 GB'

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
    --outdir ./
    """
}


workflow {
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    fastqc(ch_reads)

    if (params.run_fastp) {
        fastp(ch_reads)
        fastq_screen(fastp.out.trim_reads)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.json).collect() )

    } else {
        fastq_screen(ch_reads)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips).collect() )
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}