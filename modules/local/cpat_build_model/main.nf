process CPAT_BUILD_MODEL {
    tag "build_cpat_model"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cpat:3.0.5--py39hff726c5_3' :
        'biocontainers/cpat:3.0.5--py39hff726c5_4' }"

    input:
    path coding_fasta           // From ch_cds_ref
    path noncoding_fasta        // From ch_lncrna_ref

    output:
    path "*.Hexamer.tsv"       , emit: hexamer
    path "*.logitModel.RData"  , emit: logit_model
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "CPAT_Model"
    """
    echo "Building CPAT models from training data..."

    # Build hexamer table
    make_hexamer_tab.py \\
        -c $coding_fasta \\
        -n $noncoding_fasta \\
        > ${prefix}.Hexamer.tsv

    # Build logit model
    make_logitModel.py \\
        -c $coding_fasta \\
        -n $noncoding_fasta \\
        -x ${prefix}.Hexamer.tsv \\
        -o ${prefix}.logitModel

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpat: \$(cpat.py --version 2>&1 | grep -oP 'CPAT-\\K[0-9.]+' || echo "3.0.5")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "CPAT_Model"
    """
    touch ${prefix}.Hexamer.tsv
    touch ${prefix}.logitModel.RData

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpat: "3.0.5"
    END_VERSIONS
    """
}
