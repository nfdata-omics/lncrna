process GENERATE_LNCRNA_REPORT {
    tag "final_report"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path de_results
    path de_plots
    path classification
    path classification_stats
    path cpat_comparison
    path count_summary
    path final_gtf
    path rename_mapping

    output:
    path "*.final_report.html" , emit: report
    path "*.summary.txt"       , emit: summary
    path "*.stats.json"        , emit: stats_json
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'generate_lncrna_report.py'
}
