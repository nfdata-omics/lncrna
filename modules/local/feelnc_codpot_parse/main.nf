process FEELNC_CODPOT_PARSE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'biocontainers/biopython:1.65' }"

    input:
    tuple val(meta), path(feelnc_rf_txt)      // *_RF.txt codpot_full
    tuple val(meta), path(candidate_fasta)    // exon_filtered.fa

    output:
    tuple val(meta), path("*.feelnc.tsv") , emit: feelnc_results
    tuple val(meta), path("*.lncRNA.fa")  , emit: lncrna_fasta
    tuple val(meta), path("*.mRNA.fa")    , emit: mrna_fasta
    path "versions.yml"                   , emit: versions

    script:
    template 'parse_feelnc_output.py'
}
