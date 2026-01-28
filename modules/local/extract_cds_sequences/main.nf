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
    # Filter GTF for protein_coding transcripts
    awk -F'\\t' '(\$3=="transcript") && (\$9 ~ /(transcript_biotype|gene_biotype|gene_type) "protein_coding"/) { if (match(\$9, /transcript_id "([^"]+)"/, m)) print m[1] }' $gtf | sort -u > cds_ids.txt

    if [ -s cds_ids.txt ]; then
        # Extract entries for these transcripts
        awk 'FNR==NR{ids[\$1]=1; next} { if (match(\$0, /transcript_id "([^"]+)"/, m) && ids[m[1]]) print }' cds_ids.txt $gtf > cds.gtf

        # Extract only CDS sequences from protein-coding genes
        gffread -x cds_sequences.fa -g $fasta cds.gtf
    else
        touch cds_sequences.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
