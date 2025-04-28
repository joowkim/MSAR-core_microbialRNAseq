nextflow.enable.dsl=2

process fastp {
    tag "${meta.sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/fastp/", mode : "copy"

    module 'fastp/0.21.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.sample_name), path("${meta.sample_name}_trimmed_R{1,2}.fastq.gz"), val(meta.single_end), emit: trim_reads
    path("${meta.sample_name}.fastp.json"), emit: fastp_json

    // --trim_front2 1
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

    module "python/3.11.1"

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

    publishDir "${launchDir}/analysis/fastq_screen", mode : "copy"

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


process bowtie2 {
    tag "bowtie2 on ${sample_name}"
    label "process_high"

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
    index = params.bowtie2.(params.genome)
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
        samtools fastq -@ 6 ${sample_name}.sorted.bam \
            -1 ${sample_name}.host_remove.R1.fastq.gz \
            -2 ${sample_name}.host_remove.R2.fastq.gz \
            -0 /dev/null -s /dev/null -n
    """
}


process metaphlan {
    tag "metaphlan on ${sample_name}"
    label "memory_medium"

    publishDir "${launchDir}/analysis/metaphlan", mode : "copy"

    module "MetaPhlAn/4.0"
    // singularity

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("*")

    script:
    def bowtie2db = "/mnt/beegfs/kimj32/reference/metaphlan4/metaphlan_databases/"
    // def bowtie2db = "/mnt/beegfs/root/MetaPhlAn/" // this is from the HPC - directory is not writable
    """
        metaphlan ${reads[0]},${reads[1]} \
        --nproc ${task.cpus} \
        --input_type fastq \
        -x mpa_vOct22_CHOCOPhlAnSGB_202212 \
        --bowtie2db  ${bowtie2db} \
        -t rel_ab_w_read_stats \
        -o ${sample_name}.profiled_metagenome.txt \
        --bowtie2out ${sample_name}.bowtie2.bz2
    """
}


process concat_fq {
    // See https://forum.biobakery.org/t/humann3-paired-end-reads/862
    tag "concat ${sample_name}"
    label "process_dual"

    publishDir "${launchDir}/analysis/concat"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.concat.fastq.gz"), emit: concat_reads

    script:
    """
        cat ${reads[0]} ${reads[1]} >   ${sample_name}.concat.fastq.gz
    """
}


process megahit {
    tag "${sample_name}"
    label "memory_medium"

    publishDir (
        path: "${launchDir}/analysis/megahit",
        pattern: "**",
        saveAS: { filename -> filename.endsWith("final.contigs.fa") ? "${sample_name}/${sample_name}.contigs.fa" : filename}
    )

    module "megahit/1.2.9"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("${sample_name}", type: "dir")
    tuple val(sample_name), path("${sample_name}/final.contigs.fa"), emit: final_contigs_fa
    path("${sample_name}.*.log"), emit: "log"

    script:
    """
        megahit -1 ${reads[0]} -2 ${reads[1]} \
            -o ${sample_name} \
            -t ${task.cpus} \
            --k-min 27 \
            --k-max 47 \
            1> ${sample_name}.stdout.log \
            2> ${sample_name}.stderr.log
    """
}


process humann {
    tag "humann on ${sample_name}"
    label "memory_high"

    publishDir "${launchDir}/analysis/humann", mode : "copy"

    module "HUMAnN/3.0"
    // singularity

    input:
    tuple val(sample_name), path(reads)

    output:
    path("*")

    // --nucleotide-database : /mnt/beegfs/kimj32/reference/Humann/chocophlan
    // --protein-database : /mnt/beegfs/kimj32/reference/Humann/uniref
    // /cm/shared/apps/MetaPhlAn/4.0/envs/mpa/bin/metaphlan
    script:
    def nucl_db = "/mnt/beegfs/kimj32/reference/Humann/chocophlan"
    def prot_db = "/mnt/beegfs/kimj32/reference/Humann/uniref"
    def bowtie2_db = "/cm/shared/apps/HUMAnN/3.0/lib/python3.9/site-packages/metaphlan/metaphlan_databases"
    // def cpus = task.cpus - 2
    """
        humann --input ${reads} \
        --threads ${task.cpus} \
        --memory-use maximum \
        -o ${sample_name}_humann_out \
        --metaphlan-options "--bowtie2db ${bowtie2_db} --offline"
    """
}


process kraken2{
    tag "kraken2 on ${sample_name}"
    label "process_high"

    publishDir "${launchDir}/analysis/kraken2"

    module "kraken/2.1.2"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("*"), emit: kraken2_output

    script:
    def kraken2_db = "/mnt/beegfs/kimj32/reference/KRAKEN_DB"
    """
        kraken2 --db ${kraken2_db} \
        --threads ${task.cpus} \
        --use-names \
        --output ${sample_name}_kraken.txt \
        --report ${sample_name}_report.txt \
        --paired \
        ${reads}
     """
}


// https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
process StrainPhlAn{
    tag "StrainPhlAn on ${sample_name}"
    label "memory_medium"

    publishDir "${launchDir}/analysis/StrainPhlAn"

    module "MetaPhlAn/4.1.1"
    module 'bowtie2/2.3.4.1'
    // singularity

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    // path("*")
    tuple val(sample_name), path("${sample_name}.sam.bz2"), emit: StrainPhlAn_sam_bz2

    script:
    def bowtie2db = "/mnt/beegfs/kimj32/reference/metaphlan4/metaphlan_databases/"
    // def bowtie2db = "/mnt/beegfs/root/MetaPhlAn/" // this is from the HPC - directory is not writable
    """
        metaphlan ${reads[0]},${reads[1]} \
        --nproc ${task.cpus} \
        --input_type fastq \
        -s ${sample_name}.sam.bz2 \
        -o ${sample_name}.profiled.tsv \
        --bowtie2out ${sample_name}.bowtie2.bz2
    """
}


process concensus_markers {
    tag "concensus_markers on ${sample_name}"
    label "memory_medium"

    publishDir "${launchDir}/analysis/concensus_markers", mode : "copy"

    module "MetaPhlAn/4.0"
    module 'bowtie2/2.3.4.1'

    input:
    tuple val(sample_name), path(StrainPhlAn_sam_bz2)

    output:
    path("consensus_markers")

    script:
    """
        sample2markers.py -i ${StrainPhlAn_sam_bz2} -o consensus_markers -n ${task.cpus}
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
    // samplesheet applied // ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    fastqc(ch_reads)
    fastp(ch_reads)
    fastq_screen(fastp.out.trim_reads)
    bowtie2(fastp.out.trim_reads)
    split_reads_from_unmapped(bowtie2.out.bowtie2_bam_both_unmapped_bam)


    if (params.run_metaphlan) {
        //metaphlan(fastp.out.trim_reads)
        metaphlan(split_reads_from_unmapped.out.split_reads)
    }
    if (params.run_megahit) {
        //megahit(fastp.out.trim_reads)
        megahit(split_reads_from_unmapped.out.split_reads)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.fastp_json, bowtie2.out.samtools_stats, megahit.out.log).collect() )
    }
    if (params.run_humann) {
        concat_fq(split_reads_from_unmapped.out.split_reads)
        humann(concat_fq.out.concat_reads)
    }

    if (params.run_kraken2) {
        kraken2(split_reads_from_unmapped.out.split_reads)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.fastp_json, bowtie2.out.samtools_stats, kraken2.out.kraken2_output).collect() )
    }

    if (params.run_StrainPhlAn) {
        StrainPhlAn(split_reads_from_unmapped.out.split_reads)
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}