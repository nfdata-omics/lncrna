process FILTER_TRANSCRIPTS {
    tag "filter"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::biopython=1.79 bioconda::gffutils=0.11.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4c80ba6f0c8b13e8cf4f0c1ede8f8b6ba9d0c00d:d0e0e94bef6cdb42a25f4a5c8b42d8e7c5b5c52c-0' :
        'quay.io/biocontainers/mulled-v2-4c80ba6f0c8b13e8cf4f0c1ede8f8b6ba9d0c00d:d0e0e94bef6cdb42a25f4a5c8b42d8e7c5b5c52c-0' }"

    input:
    path gtf
    val min_length
    val min_exons

    output:
    path "filtered_transcripts.gtf"   , emit: gtf
    path "filtered_transcripts.fasta" , emit: fasta
    path "filter_stats.txt"           , emit: stats
    path "versions.yml"               , emit: versions

    script:
    template 'filter_transcripts.py'
}
