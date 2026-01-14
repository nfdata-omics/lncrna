process EXTRACT_MRNA_SEQUENCES {
    tag "extract_mrna"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"
    input:
    path fasta
    path gtf

    output:
    path "mrna_sequences.fa", emit: fasta
    path "versions.yml"     , emit: versions

    script:
    """
    # Extract full mRNA transcripts from protein-coding genes
    gffread -w mrna_sequences.fa -g $fasta $gtf -C

    # Filter to keep only protein_coding transcripts
    grep -A 1 'transcript_biotype "protein_coding"' mrna_sequences.fa > mrna_sequences_filtered.fa || true

    # If filtering didn't work (no biotype), use all transcripts
    if [ ! -s mrna_sequences_filtered.fa ]; then
        mv mrna_sequences.fa mrna_sequences_filtered.fa
    fi
    mv mrna_sequences_filtered.fa mrna_sequences.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
