process FILTER_TRANSCRIPTS_LENGTH {
    tag "length_filter"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(gtf)                      // From ch_gtf
    tuple val(meta2), path(fasta)                   // From ch_transcripts_fa

    output:
    tuple val(meta), path("*.length_filtered.gtf")  , emit: filtered_length_gtf
    tuple val(meta2), path("*.length_filtered.fa")  , emit: filtered_length_fasta
    path "*.length_stats.txt"                       , emit: stats
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_length = task.ext.args ?: '200'
    """
    filter_transcripts_by_length.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        --min_length $min_length \\
        --prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
