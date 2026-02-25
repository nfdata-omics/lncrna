process GTF_FILTER {
    tag "$fasta"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path fasta
    path gtf

    output:
    path "*.filtered.gtf", emit: genome_gtf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'filter_gtf.py'

    stub:
    """
    touch ${fasta.baseName}.filtered.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
