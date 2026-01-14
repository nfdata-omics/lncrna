process RENAME_LNCRNA {
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
    tuple val(meta), path("*.renamed.gtf"), emit: gtf
    path "*.mapping.txt"                   , emit: mapping
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rename_lncrna.py \\
        --lncrna_gtf $lncrna_gtf \\
        --protein_gtf $protein_coding_gtf \\
        --output ${prefix}.renamed.gtf \\
        --mapping ${prefix}.mapping.txt \\
        --prefix LINC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
