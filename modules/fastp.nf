process fastp {
    tag "${meta.sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/fastp/"

    module "fastp/1.0.1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.sample_name), path("${meta.sample_name}_trimmed_R{1,2}.fastq.gz"), val(meta.single_end), emit: trim_reads
    path("${meta.sample_name}.fastp.json"), emit: fastp_json

    script:
    if (!meta.single_end) {
        """
        fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        --thread ${task.cpus} \
        -o ${meta.sample_name}_trimmed_R1.fastq.gz \
        -O ${meta.sample_name}_trimmed_R2.fastq.gz \
        --json ${meta.sample_name}.fastp.json
        """
    } else {
        """
        fastp \
        -i ${reads} \
        --thread ${task.cpus} \
        -o ${meta.sample_name}_trimmed_R1.fastq.gz \
        --json ${meta.sample_name}.fastp.json
        """
    }
}
