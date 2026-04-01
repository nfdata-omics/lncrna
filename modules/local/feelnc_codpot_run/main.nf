process FEELNC_CODPOT_RUN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/feelnc:0.2--pl526_0' :
        'biocontainers/feelnc:0.2--pl526_0' }"

    input:
    tuple val(meta), path(candidate_fasta)      // FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
    path mrna_fasta                             // PREPARE_GENOME.out.mrna_fasta

    output:
    tuple val(meta), path("*.feelnc_codpot_RF.txt") , emit: codpot_full
    tuple val(meta), path("*.feelnc_codpot.noORF.fa"), emit: candidate_fasta
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ""
    """
    export FEELNCPATH=\$(dirname \$(which FEELnc_codpot.pl))/..

    FEELnc_codpot.pl \\
        -i ${candidate_fasta} \\
        -a ${mrna_fasta} \\
        --mode=shuffle \\
        --numtx=500,500 \\
        -o ${prefix}.feelnc_codpot \\
        ${args}

    ls -la feelnc_codpot_out/ || true

    mv feelnc_codpot_out/${prefix}.feelnc_codpot_RF.txt .
    mv feelnc_codpot_out/${prefix}.feelnc_codpot.noORF.fa .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        feelnc: "0.2.1"
    END_VERSIONS
    """
}
