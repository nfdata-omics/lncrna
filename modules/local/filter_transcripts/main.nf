process FILTER_TRANSCRIPTS {
    tag "filter"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::biopython=1.79 bioconda::gffutils=0.11.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4c80ba6f0c8b13e8cf4f0c1ede8f8b6ba9d0c00d:d0e0e94bef6cdb42a25f4a5c8b42d8e7c5b5c52c-0' :
        'quay.io/biocontainers/mulled-v2-4c80ba6f0c8b13e8cf4f0c1ede8f8b6ba9d0c00d:d0e0e94bef6cdb42a25f4a5c8b42d8e7c5b5c52c-0' }"

    input:
    path gtf
    val min_length
    val min_exons

    output:
    path "filtered_transcripts.gtf"   , emit: gtf
    path "filtered_transcripts.fasta" , emit: fasta
    path "filter_stats.txt"           , emit: stats
    path "versions.yml"               , emit: versions

    script:
    """
    #!/usr/bin/env python3

    import gffutils
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import sys

    # Create database from GTF
    print("Creating GTF database...")
    db = gffutils.create_db('$gtf', 'transcripts.db', force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)

    db = gffutils.FeatureDB('transcripts.db')

    # Filter transcripts
    filtered_count = 0
    total_count = 0
    filtered_transcripts = []

    print("Filtering transcripts...")
    with open('filtered_transcripts.gtf', 'w') as out_gtf, \\
        open('filtered_transcripts.fasta', 'w') as out_fasta:

        for transcript in db.features_of_type('transcript'):
            total_count += 1

            # Get exons
            exons = list(db.children(transcript, featuretype='exon'))
            num_exons = len(exons)

            # Calculate transcript length
            transcript_length = sum([exon.end - exon.start + 1 for exon in exons])

            # Apply filters
            if transcript_length >= ${min_length} and num_exons >= ${min_exons}:
                filtered_count += 1

                # Write GTF
                out_gtf.write(str(transcript) + '\\n')
                for exon in exons:
                    out_gtf.write(str(exon) + '\\n')

                # Create FASTA entry (placeholder sequence)
                # In real analysis, extract from genome
                transcript_id = transcript.attributes.get('transcript_id', [transcript.id])[0]
                seq_record = SeqRecord(
                    Seq('N' * transcript_length),
                    id=transcript_id,
                    description=f"length={transcript_length} exons={num_exons}"
                )
                SeqIO.write(seq_record, out_fasta, 'fasta')

    # Write statistics
    with open('filter_stats.txt', 'w') as stats:
        stats.write(f"Total transcripts: {total_count}\\n")
        stats.write(f"Filtered transcripts: {filtered_count}\\n")
        stats.write(f"Percentage retained: {100*filtered_count/total_count:.2f}%\\n")
        stats.write(f"Min length: ${min_length}\\n")
        stats.write(f"Min exons: ${min_exons}\\n")

    print(f"Filtered {filtered_count} out of {total_count} transcripts")

    # Write versions
    with open('versions.yml', 'w') as v:
        v.write('"${task.process}":\\n')
        v.write('    python: "' + sys.version.split()[0] + '"\\n')
    """
}
