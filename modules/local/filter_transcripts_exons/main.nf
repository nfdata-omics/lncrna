process FILTER_TRANSCRIPTS_EXONS {
    tag "exon_filter"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf)                  // From FILTER_TRANSCRIPTS_LENGTH.out.filtered_length_gtf
    tuple val(meta2), path(fasta)               // From FILTER_TRANSCRIPTS_LENGTH.out.filtered_length_fasta

    output:
    tuple val(meta), path("*.exon_filtered.gtf")  , emit: filtered_exon_gtf
    tuple val(meta2), path("*.exon_filtered.fa")  , emit: filtered_exon_fasta
    path "*.exon_stats.txt"                       , emit: stats
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_exons = task.ext.args ?: '2'
    """
    filter_transcripts_by_exons.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        --min_exons $min_exons \\
        --prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
