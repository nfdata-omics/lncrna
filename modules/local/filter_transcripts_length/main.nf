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
    template 'filter_transcripts_by_length.py'
}
