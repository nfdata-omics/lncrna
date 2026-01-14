process EXTRACT_LNCRNA_SEQUENCES {
    tag "extract_lncrna"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

    input:
    path fasta
    path gtf

    output:
    path "lncrna_sequences.fa", emit: fasta
    path "versions.yml"        , emit: versions

    script:
    """
    # Extract lncRNA transcripts
    gffread -w all_transcripts.fa -g $fasta $gtf

    # Filter to keep only lncRNA biotypes
    grep -A 1 -E 'transcript_biotype "(lncRNA|lincRNA|antisense|sense_intronic|sense_overlapping)"' \\
        all_transcripts.fa > lncrna_sequences.fa || true

    # If no lncRNAs found, create empty file
    if [ ! -s lncrna_sequences.fa ]; then
        touch lncrna_sequences.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
