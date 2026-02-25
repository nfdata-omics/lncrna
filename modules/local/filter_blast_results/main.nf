process FILTER_BLAST_RESULTS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4' :
        'biocontainers/blast:2.17.0--h66d330f_0' }"

    input:
    tuple val(meta), path(lncrna_fasta)         // tuple val(meta), path(fasta)
    tuple val(meta), path(blast_hits)           // tuple val(meta), path(blast_hits)

    output:
    tuple val(meta), path("*.filtered.fa")     , emit: filtered_fasta
    path "*.blast_filter_stats.txt"            , emit: stats
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'filter_blast_results.py'
}
