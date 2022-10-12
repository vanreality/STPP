process BUILD_INTERVALS {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta_fai)

    output:
    tuple val(meta), path("${fasta_fai.baseName}.bed") , emit: bed
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fasta_fai} > ${fasta_fai.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
