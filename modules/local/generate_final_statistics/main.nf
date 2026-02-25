process GENERATE_FINAL_STATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(final_gtf)
    tuple val(meta2), path(lncrna_gtf)
    path cpat_lncrna_results
    path cpat_coding_results
    path classification

    output:
    path "*.final_report.html"          , emit: final_stats_report
    path "*.statistics_summary.tsv"     , emit: final_stats_summary
    path "*.cpat_comparison.png"        , emit: cpat_plot
    path "*.classification_barplot.png" , emit: cpat_classification_plot
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'generate_final_statistics.py'
}
