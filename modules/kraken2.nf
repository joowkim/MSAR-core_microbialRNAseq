process kraken2 {
    tag "kraken2 on ${sample_name}"
    label "process_high"

    publishDir "${launchDir}/analysis/kraken2", mode: "copy"

    module "kraken/2.1.2"

    input:
    tuple val(sample_name), path(reads)

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
