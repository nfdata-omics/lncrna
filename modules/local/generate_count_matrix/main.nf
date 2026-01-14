process GENERATE_COUNT_MATRIX {
    tag "count_matrix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path count_files  // All count files from samples
    path gtf          // GTF annotation
    val method        // Quantification method used

    output:
    path "*.count_matrix.tsv"    , emit: matrix
    path "*.gene_info.tsv"       , emit: gene_info
    path "*.sample_summary.tsv"  , emit: summary
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "expression"
    """
    merge_count_matrix.py \\
        --counts ${count_files} \\
        --gtf $gtf \\
        --method $method \\
        --output ${prefix}.count_matrix.tsv \\
        --gene_info ${prefix}.gene_info.tsv \\
        --summary ${prefix}.sample_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
