process kraken_biom {
    tag "kraken_biom"
    label "process_low"

    publishDir "${launchDir}/analysis/kraken_biom", mode: "copy"

    input:
    path(reports)

    output:
    path("biom_table.tsv"), emit: biom_table

    script:
    def kraken_biom = "~/beegfs/python_env/kraken2/bin/kraken-biom"
    """
    ${kraken_biom} ${reports} --fmt tsv -o biom_table.tsv
    """
}
