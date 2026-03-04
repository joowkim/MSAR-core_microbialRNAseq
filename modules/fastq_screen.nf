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

    // threads option is already defined in fastq_screen.conf
    script:
    def conf = "/mnt/beegfs/kimj32/polymerase/polymeraseDependencies/FastQ_Screen_Genomes/fastq_screen.conf"
    """
    fastq_screen --aligner bowtie2 \
    --conf ${conf} \
    ${reads[0]} \
    --outdir ./ \
    --threads ${task.cpus}
    """
}
