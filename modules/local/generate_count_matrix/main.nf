process GENERATE_COUNT_MATRIX {
    tag "count_matrix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path count_files  // All count files from samples
    path gtf          // GTF annotation
    val method        // Quantification method used
    val lncrna_biotypes

    output:
    path "*.count_matrix.tsv"    , emit: matrix
    path "*.lncrna_matrix.tsv"   , emit: lncrna_matrix
    path "*.gene_info.tsv"       , emit: gene_info
    path "*.sample_summary.tsv"  , emit: summary
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'merge_count_matrix.py'
}
