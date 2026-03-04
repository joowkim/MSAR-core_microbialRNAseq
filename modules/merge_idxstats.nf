process merge_idxstats {
    tag "Merging idxstats"
    label "process_low"

    publishDir "${launchDir}/analysis/gene_count_matrix", mode: "copy"

    input:
    path(idxstats)

    output:
    path("gene_count_matrix.txt")

    script:
    """
    first=\$(ls *.idxstat | sort | head -n 1)
    cut -f1 "\$first" > genes.txt

    for f in \$(ls *.idxstat | sort); do
        cut -f3 "\$f" > "\${f}.count"
    done

    echo -e "Gene\\t\$(ls *.idxstat | sort | sed 's/.idxstat//g' | tr '\\n' '\\t')" > header.txt
    paste genes.txt \$(ls *.idxstat | sort | sed 's/\$/.count/') > body.txt
    cat header.txt body.txt > gene_count_matrix.txt
    """
}
