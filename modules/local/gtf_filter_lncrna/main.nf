process GTF_FILTER_LNCRNA {
    tag "$gtf"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path gtf
    val lncrna_biotypes

    output:
    path "*.lncrna.gtf"   , emit: lncrna_gtf
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'filter_gtf_lncrna.py'

    stub:
    """
    touch ${gtf.baseName}.lncrna.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
