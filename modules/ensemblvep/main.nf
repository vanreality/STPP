process ENSEMBLVEP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf)
    val   genome
    val   species
    val   cache_version
    val   cache
    path  fasta
    path  extra_files

    output:
    tuple val(meta), path("*.ann.vcf")  , optional:true, emit: vcf
    tuple val(meta), path("*.ann.tab")  , optional:true, emit: tab
    tuple val(meta), path("*.ann.json") , optional:true, emit: json
    path "*.summary.html"               , emit: report
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""

    """
    vep \\
        -i $vcf \\
        -o ${prefix}.ann.${file_extension} \\
        $args \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
