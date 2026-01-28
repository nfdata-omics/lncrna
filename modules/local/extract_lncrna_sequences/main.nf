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
    awk -F'\\t' '(${'$'}3=="transcript") && (${'$'}9 ~ /(transcript_biotype|gene_biotype|gene_type) "(lncRNA|lincRNA|antisense|sense_intronic|sense_overlapping|processed_transcript|non_coding|3prime_overlapping_ncRNA|non_coding_gene)"/) { if (match(${'$'}9, /transcript_id "([^"]+)"/, m)) print m[1] }' $gtf | sort -u > lnc_ids.txt
    if [ -s lnc_ids.txt ]; then
        # Extract entries for these transcripts, excluding CDS/codon features to avoid CDS= in headers
        awk 'FNR==NR{ids[${'$'}1]=1; next} { if (match(${'$'}0, /transcript_id "([^"]+)"/, m) && ids[m[1]] && ${'$'}3!="CDS" && ${'$'}3!="start_codon" && ${'$'}3!="stop_codon") print }' lnc_ids.txt $gtf > lnc.gtf
    else
        > lnc.gtf
    fi
    if [ -s lnc.gtf ]; then
        gffread -w lncrna_sequences.fa -g $fasta lnc.gtf
    else
        > lncrna_sequences.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
