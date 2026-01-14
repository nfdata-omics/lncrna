process EXTRACT_CDS_SEQUENCES {
    tag "extract_cds"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

    input:
    path fasta
    path gtf

    output:
    path "cds_sequences.fa", emit: fasta
    path "versions.yml"    , emit: versions

    script:
    """
    # Extract only CDS sequences from protein-coding genes
    gffread -x cds_sequences.fa -g $fasta $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
