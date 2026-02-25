process FEELNC_CODPOT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/feelnc:0.2--pl526_0' :
        'biocontainers/feelnc:0.2--pl526_0' }"

    input:
    tuple val(meta), path(candidate_fasta)        // From ch_filtered_exons_fa
    path genome_fasta                             // From ch_fasta

    output:
    tuple val(meta), path("*.feelnc.tsv")          , emit: feelnc_results
    path "*.feelnc_codpot.txt"                     , emit: codpot_full
    path "*.lncRNA.fa"                             , emit: lncrna_fasta, optional: true
    path "*.mRNA.fa"                               , emit: mrna_fasta, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'parse_feelnc_output.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.feelnc.tsv
    touch ${prefix}.feelnc_codpot.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        feelnc: "0.2.1"
    END_VERSIONS
    """
}
