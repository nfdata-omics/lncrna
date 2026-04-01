process PLEK_PARSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plek:1.2--py39he88f293_9' :
        'biocontainers/plek:1.2--py311h8ddd9a4_10' }"

    input:
    //tuple val(meta), path(fasta)
    tuple val(meta), path(plek_raw)

    output:
    tuple val(meta), path("*.plek.tsv"), emit: plek_results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'convert_plek_output.py'
}
