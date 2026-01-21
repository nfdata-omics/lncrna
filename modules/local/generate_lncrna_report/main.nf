process GENERATE_LNCRNA_REPORT {
    tag "final_report"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

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
    def prefix = task.ext.prefix ?: "lncRNA_pipeline"
    def de_results_arg = de_results.name != 'NO_FILE' ? "--de_results $de_results" : ""
    def de_plots_arg = de_plots.name != 'NO_FILE' ? "--de_plots $de_plots" : ""
    def count_summary_arg = count_summary.name != 'NO_FILE' ? "--count_summary $count_summary" : ""
    def rename_mapping_arg = rename_mapping.name != 'NO_FILE' ? "--rename_mapping $rename_mapping" : ""
    def cpat_arg = cpat_comparison.name != 'NO_FILE' ? "--cpat_comparison $cpat_comparison" : "--cpat_comparison NO_FILE"
    """
    generate_lncrna_report.py \\
        $de_results_arg \\
        $de_plots_arg \\
        --classification $classification \\
        --classification_stats $classification_stats \\
        $cpat_arg \\
        $count_summary_arg \\
        --final_gtf $final_gtf \\
        $rename_mapping_arg \\
        --output ${prefix}.final_report.html \\
        --summary ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "not available")
    END_VERSIONS
    """
}
