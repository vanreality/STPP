process GATK4_MERGEVCFS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf)
    path  dict

    output:
    tuple val(meta), path('*.vcf.gz'), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = vcf.collect{ "--INPUT $it"}.join(' ')
    def reference_command = dict ? "--SEQUENCE_DICTIONARY $dict" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MergeVcfs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" MergeVcfs \\
        $input_list \\
        --OUTPUT ${prefix}.vcf.gz \\
        $reference_command \\
        --TMP_DIR . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
