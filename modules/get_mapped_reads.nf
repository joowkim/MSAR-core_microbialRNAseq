process get_mapped_reads {
    tag "get_mapped_reads ${sample_name}"
    label "process_low"

    publishDir "${launchDir}/analysis/mapped_reads_from_bam"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam_file)

    output:
    tuple val(sample_name), path("${sample_name}_R{1,2}.fastq.gz"), emit: mapped_reads

    script:
    """
    samtools fastq ${bam_file} \
        -1 ${sample_name}_R1.fastq.gz \
        -2 ${sample_name}_R2.fastq.gz \
        -0 /dev/null -s /dev/null -n
    """
}
