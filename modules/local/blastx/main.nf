process BLASTX {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4' :
        'biocontainers/blast:2.17.0--h66d330f_0' }"

    input:
    tuple val(meta), path(lncrna_fasta)
    path blast_protein_db       // BLAST database of protein sequences
    val evalue_cutoff

    output:
    tuple val(meta), path("*.blast_hits.txt")       , emit: blast_hits
    tuple val(meta), path("*.removed_ids.txt")      , emit: removed_ids
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def evalue = evalue_cutoff ?: '1e-5'
    """
    # Run BLASTx against protein database
    blastx \\
        -query $lncrna_fasta \\
        -db $blast_protein_db \\
        -evalue $evalue \\
        -num_threads $task.cpus \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
        -max_target_seqs 1 \\
        $args \\
        -out ${prefix}.blast_hits.txt

    # Extract IDs with significant hits (to remove)
    if [ -s ${prefix}.blast_hits.txt ]; then
        cut -f1 ${prefix}.blast_hits.txt | sort -u > ${prefix}.removed_ids.txt
    else
        touch ${prefix}.removed_ids.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
