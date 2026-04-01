process DIFFERENTIAL_EXPRESSION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker.io/nfdata/edger:v4.4.2"

    input:
    tuple val(meta), path(count_matrix)
    path design_file
    val   contrast
    val   LFC
    val   FDR
    val   min_count
    val   min_total
    val   correction_mode
    val   batch
    val   ruv_k
    val   ruv_controls
    val   top_genes
    val de_method

    output:
    tuple val(meta), path("DE_*/*DGE_QLF_*.xlsx")          , emit: results_xlsx, optional: true
    tuple val(meta), path("DE_*/*DGE_LRT_*.txt")           , emit: results_lrt , optional: true
    tuple val(meta), path("DE_*/*normalized_CPM_*.csv")    , emit: normalized  , optional: true
    tuple val(meta), path("DE_*/*DGE_summary_*.txt")       , emit: summary     , optional: true
    tuple val(meta), path("DE_*/*DGE_summary_*.tsv")       , emit: summary_tsv     , optional: true
    tuple val(meta), path("DE_*/plots/png/*.png")          , emit: plots_png   , optional: true
    tuple val(meta), path("DE_*/plots/pdf/*.pdf")          , emit: plots_pdf   , optional: true
    tuple val(meta), path("DE_*/MDS_mqc.png")              , emit: mds_png     , optional: true
    tuple val(meta), path("DE_*/PCA_mqc.png")              , emit: pca_png     , optional: true
    tuple val(meta), path("DE_*/heatmap_global_mqc.png")   , emit: heatmap_png , optional: true
    tuple val(meta), path("DE_*/MDS.pdf")                  , emit: mds_pdf     , optional: true
    tuple val(meta), path("DE_*/PCA.pdf")                  , emit: pca_pdf     , optional: true
    tuple val(meta), path("DE_*/heatmap_global.pdf")       , emit: heatmap_pdf , optional: true
    tuple val(meta), path("DE_*/DGE_lncRNAs.RData")        , emit: rdata       , optional: true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def contrast        = task.ext.contrast ?: (meta.contrast ?: "")
    def LFC             = task.ext.LFC ?: ""
    def FDR             = task.ext.FDR ?: ""
    def min_count        = task.ext.min_count ?: ""
    def min_total        = task.ext.min_total ?: ""
    def correction_mode  = task.ext.correction_mode ?: ""
    def batch            = task.ext.batch ?: ""
    def ruv_k            = task.ext.ruv_k ?: ""
    def ruv_controls     = task.ext.ruv_controls ?: ""
    def top_genes        = task.ext.top_genes ?: ""

    template 'differential_expression.R'

    stub:
    """
    mkdir -p DE_stub/plots/png DE_stub/plots/pdf
    touch DE_stub/DGE_QLF_stub.xlsx
    touch DE_stub/DGE_LRT_stub.txt
    touch DE_stub/normalized_CPM_stub.csv
    touch DE_stub/DGE_summary_stub.txt
    touch DE_stub/plots/png/stub.png
    touch DE_stub/plots/pdf/stub.pdf

    touch MDS.png PCA.png heatmap_global.png
    touch MDS.pdf PCA.pdf heatmap_global.pdf
    touch DGE_lncRNAs.RData
    """
}
