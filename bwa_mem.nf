nextflow.enable.dsl=2


process fastp {
    tag "${meta.sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/fastp/"

    module 'fastp/0.21.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.sample_name), path("${meta.sample_name}_trimmed_R{1,2}.fastq.gz"), val(meta.single_end), emit: trim_reads
    path("${meta.sample_name}.fastp.json"), emit: fastp_json

    script:
    def adapter = "/mnt/beegfs/kimj32/reference/adapters.fa"
    if(!meta.single_end) {
    """
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    --thread ${task.cpus} \
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
    -o ${meta.sample_name}_trimmed_R1.fastq.gz \
    --adapter_fasta $adapter \
    --json ${meta.sample_name}.fastp.json
    """
    }
}


process fastqc {
    tag "Fastqc on ${meta.sample_name}"
    label "process_low"

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
    tag "Multiqc on this project"
    label "process_dual"

    publishDir "${launchDir}/analysis/multiqc/", mode : "copy"

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    config_yaml = "/home/kimj32/config_defaults.yaml"
    """
    multiqc ${files} --filename "multiqc_report.html" --config ${config_yaml}
    """
}


process fastq_screen {
    tag "Fastq-screen on ${sample_name}"
    label "process_low"

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
    tag "${sample_name}"
   label "process_medium"

    publishDir "${launchDir}/analysis/bwa_mem"

    module "bwa/0.7.17"
    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), emit: bam_file

    // -M: mark shorter split hits as secondary
    script:
    def index = params.bwa_index
    //def index = "/mnt/beegfs/kimj32/reference/Marker_genes/cutC_D.fasta"
    // def index = "/mnt/beegfs/kimj32/reference/Marker_genes/gbuA.fasta"
    """
    bwa mem -M \
    -t ${task.cpus} \
     ${index} \
     ${reads} \
     | samtools sort \
     -O "BAM" \
     -o ${sample_name}.bam -
    """
}

process samtools {
    tag "${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/samtools_stat", mode:"copy"
    module "samtools/1.16.1"

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

process bowtie2 {
    tag "bowtie2 on ${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/bowtie2"

    module 'bowtie2/2.3.4.1'
    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.mapped_unmapped.bam"), emit: bowtie2_mapped_unmapped_bam
    tuple val(sample_name), path("${sample_name}.both_unmapped.bam"), val(is_SE), emit: bowtie2_bam_both_unmapped_bam
    path("${sample_name}.mapped_unmapped.stats"), emit: samtools_stats

    script:
    index = "/mnt/beegfs/kimj32/reference/mouse/gencode/GRCm38.p6/indexes/bowtie2/mouse"
    """
    bowtie2 -p ${task.cpus} -x ${index} \
    -1 ${reads[0]} -2 ${reads[1]} \
    | samtools sort -O BAM -o ${sample_name}.mapped_unmapped.bam

    samtools index -@ 4 ${sample_name}.mapped_unmapped.bam
    samtools stats -@ 4 ${sample_name}.mapped_unmapped.bam > ${sample_name}.mapped_unmapped.stats

    samtools view -@ 4 -b -f 12 -F 256 \
    ${sample_name}.mapped_unmapped.bam > ${sample_name}.both_unmapped.bam

    """
}


process split_reads_from_unmapped {
    tag "split reads - ${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/split_reads"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam_file), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.host_remove.R{1,2}.fastq.gz"), val(is_SE),  emit: split_reads

    script:
    """
    samtools sort -n ${bam_file} -o ${sample_name}.sorted.bam
    samtools fastq -@ 10 ${sample_name}.sorted.bam \
        -1 ${sample_name}.host_remove.R1.fastq.gz \
        -2 ${sample_name}.host_remove.R2.fastq.gz \
        -0 /dev/null -s /dev/null -n
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
    // bowtie2(fastp.out.trim_reads)
    // split_reads_from_unmapped(bowtie2.out.bowtie2_bam_both_unmapped_bam)
    // ch_reference = Channel.value(params.index)
    // bwa_mem(split_reads_from_unmapped.out.split_reads)
    bwa_mem(fastp.out.trim_reads)
    samtools(bwa_mem.out.bam_file)
    multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.fastp_json, samtools.out.stats_out).collect() )
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}