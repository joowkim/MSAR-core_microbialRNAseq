// Keep only properly paired mapped reads
process filt_bam {
    tag "filt_bam ${sample_name}"
    label "process_low"

    publishDir "${launchDir}/analysis/filt_bam"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam_file)

    output:
    tuple val(sample_name), path("${sample_name}.sorted.mapped.bam"), emit: filt_bam

    // view -f 2, keep properly paired
    script:
    """
    samtools view -b -f 2 ${bam_file} > ${sample_name}_unsorted.mapped.bam
    samtools sort -n ${sample_name}_unsorted.mapped.bam -o ${sample_name}.sorted.mapped.bam
    """
}
