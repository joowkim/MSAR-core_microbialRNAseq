nextflow.enable.dsl=2


process fastp {
    debug true
    tag "${meta.sample_name}"
    // label "universal"
    cpus 12
    memory '32 GB'

    publishDir "${launchDir}/analysis/fastp/"

    module 'fastp/0.21.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.sample_name), path("${meta.sample_name}_trimmed_R*.fastq.gz"), val(meta.single_end), emit: trim_reads
    path("${meta.sample_name}.fastp.json"), emit: fastp_json

    script:
    def adapter = "/mnt/beegfs/kimj32/reference/adapters.fa"
    if(!meta.single_end) {
    """
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    --thread ${task.cpus} \
    --qualified_quality_phred 20 \
    -o ${meta.sample_name}_trimmed_R1.fastq.gz \
    -O ${meta.sample_name}_trimmed_R2.fastq.gz \
    --adapter_fasta $adapter \
    --json ${meta.sample_name}.fastp.json
    """
    } else {
    """
    fastp \
    -i ${reads} \
    --thread ${task.cpus} \
    --qualified_quality_phred 20 \
    -o ${meta.sample_name}_trimmed_R1.fastq.gz \
    --adapter_fasta $adapter \
    --json ${meta.sample_name}.fastp.json
    """
    }
}


process fastqc {
    debug true
    tag "Fastqc on ${meta.sample_name}"
    //label "universal" // getting cpu and memory usage from salmon.config - called universal
    cpus 12
    memory '32 GB'

    publishDir "${launchDir}/analysis/fastqc/", mode: "copy"

    module 'FastQC/0.11.9'

    input:
    tuple val(meta), path(reads)

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
    multiqc ${files} --filename "multiqc_report.metaphlan4.html" --config ${config_yaml}
    """
}


process fastq_screen {
    debug true
    tag "Fastq-screen on ${sample_name}"

    cpus 8
    memory '16 GB'

    publishDir "${launchDir}/analysis/fastq_screen"

    module 'FastQScreen/0.14.1'
    module 'bowtie2/2.3.4.1'

    input:
    tuple val(sample_name), path(reads), val(is_SE)

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


process bwa_mem {
    debug true
    tag "${sample_name}"

    cpus 10
    memory '32 GB'

    publishDir "${launchDir}/analysis/bwa_mem"

    module "bwa/0.7.17"
    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), emit: bam_file

    // -M: mark shorter split hits as secondary
    script:
    def index = "/mnt/beegfs/kimj32/reference/Marker_genes/cutC_D.fasta"
    // def index = "/mnt/beegfs/kimj32/reference/Marker_genes/gbuA.fasta"
    """
    bwa mem -M \
    -t ${task.cpus} \
     ${index} \
     ${reads} \
     | samtools sort \
     -@ 4 \
     -m 6G \
     -O "BAM" \
     -o ${sample_name}.bam -
    """
}

process samtools {
    debug true
    tag "${sample_name}"

    publishDir "${launchDir}/analysis/samtools_stat", mode:"copy"

    module "samtools/1.16.1"

    cpus 12
    memory '16 GB'

    input:
    tuple val(sample_name), path(bam)

    output:
    path("${sample_name}.*stat*"), emit: stats_out
    path("${sample_name}.bam.bai")


    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools idxstats ${bam} > ${sample_name}.idxstat
    samtools stats ${bam} > ${sample_name}.stats
    samtools flagstat -O tsv ${bam} > ${sample_name}.flagstat
    """

}


// See https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

// adapted from https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_reads = ch_samplesheet.splitCsv(header:true).map {

    // This is the read1 and read2 entry
    r1 = it['fq1']
    r2 = it['fq2']

    // Detect wiether single-end or paired-end
    is_singleEnd = r2.toString() =='' ? true : false

    // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
    meta = [sample_name: it['sample'], single_end: is_singleEnd]

    // We return a nested map, the first entry is the meta map, the second one is the read(s)
    r2.toString()=='' ? [meta, [r1]] : [meta, [r1, r2]]

}


workflow {
    fastqc(ch_reads)
    fastp(ch_reads)
    fastq_screen(fastp.out.trim_reads)
    // ch_reference = Channel.value(params.index)
    bwa_mem(fastp.out.trim_reads)
    samtools(bwa_mem.out.bam_file)
    multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.fastp_json, samtools.out.stats_out).collect() )
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}