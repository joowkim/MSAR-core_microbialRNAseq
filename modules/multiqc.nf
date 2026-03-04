process multiqc {
    tag "Multiqc on this project"
    label "process_dual"

    publishDir "${launchDir}/analysis/multiqc/", mode: "copy"

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    def report_name = params.run_name ?: "multiqc_report"
    """
    multiqc ${files} --filename "${report_name}.html"
    """
}
