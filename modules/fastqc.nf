process fastqc {
    tag "Fastqc on ${meta.sample_name}"
    label "process_low"

    publishDir "${launchDir}/analysis/fastqc/", mode: "copy"

    module 'FastQC/0.11.9'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.zip"), emit: zips
    path("*.html"), emit: htmls

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
