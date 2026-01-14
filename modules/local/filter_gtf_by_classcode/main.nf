process FILTER_GTF_BY_CLASSCODE {
    tag "$gtf"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path gtf

    output:
    path "*.filtered.gtf"       , emit: filtered_gtf
    path "*.classcode_stats.txt", emit: stats
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    def classcodes = task.ext.args ?: 'i,u,x'  // Default: i, u, x
    """
    filter_gtf_by_classcode.py \\
        --gtf $gtf \\
        --classcodes $classcodes \\
        --prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.filtered.gtf
    touch ${prefix}.classcode_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
