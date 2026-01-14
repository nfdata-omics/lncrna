process GTF_FILTER_LNCRNA {
    tag "$gtf"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path gtf

    output:
    path "*.lncrna.gtf"   , emit: lncrna_gtf
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // filter_gtf_lncrna.py is our custom script - filters to lncRNA only
    """
    filter_gtf_lncrna.py \\
        --gtf $gtf \\
        --prefix ${gtf.baseName} \\
        --log_file ${gtf.baseName}.filter.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${gtf.baseName}.lncrna.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
