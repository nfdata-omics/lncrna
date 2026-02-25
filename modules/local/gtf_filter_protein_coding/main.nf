process GTF_FILTER_PROTEIN_CODING {
    tag "$gtf"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path gtf

    output:
    path "*.protein_coding.gtf" , emit: protein_gtf
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'filter_gtf_protein_coding.py'

    stub:
    """
    touch ${gtf.baseName}.protein_coding.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
