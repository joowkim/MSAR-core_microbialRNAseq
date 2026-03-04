process samtools {
    tag "${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/samtools_stat", mode: "copy"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam)

    output:
    path("${sample_name}.stats"), emit: stats_out
    path("${sample_name}.idxstat"), emit: idxstats_out
    path("${sample_name}.bam.bai")

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools idxstats ${bam} > ${sample_name}.idxstat
    samtools stats ${bam} > ${sample_name}.stats
    samtools flagstat -O tsv ${bam} > ${sample_name}.flagstat
    """
}
