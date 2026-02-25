process CLASSIFY_LNCRNA {
tag "$meta.id"
label 'process_single'

conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'biocontainers/python:3.9--1' }"

input:
tuple val(meta), path(lncrna_gtf)
path protein_coding_gtf

output:
tuple val(meta), path("*.classification.txt"), emit: classification
tuple val(meta), path("*.stats.txt")         , emit: stats
tuple val(meta), path("README.txt")          , emit: readme
path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'classify_lncrna.py'
}
