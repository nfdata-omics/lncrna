process ANALYSIS_CIS_TRANS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(expression_matrix)
    path gtf

    output:
    path "*.cis_results.tsv", emit: cis_results
    path "*.trans_results.tsv", optional: true, emit: trans_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${workflow.projectDir}/bin/analyze_cis_trans.py \\
        --expression $expression_matrix \\
        --gtf $gtf \\
        --output_cis ${prefix}.cis_results.tsv \\
        --output_trans ${prefix}.trans_results.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "not available")
        scipy: \$(python -c "import scipy; print(scipy.__version__)" 2>/dev/null || echo "not available")
    END_VERSIONS
    """
}
