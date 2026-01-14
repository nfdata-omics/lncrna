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
path "versions.yml"                          , emit: versions

when:
task.ext.when == null || task.ext.when

script:
def prefix = task.ext.prefix ?: "${meta.id}"
"""
classify_lncrna.py \\
    --lncrna_gtf $lncrna_gtf \\
    --protein_gtf $protein_coding_gtf \\
    --output ${prefix}.classification.txt \\
    --stats ${prefix}.stats.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
END_VERSIONS
"""
}
