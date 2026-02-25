#!/usr/bin/env python3
import argparse
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import platform

def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def main():
    parser = argparse.ArgumentParser(description='Filter transcripts by length and exon count, emit GTF/FASTA and stats')
    parser.add_argument('--gtf', required=True, help='Input GTF file')
    parser.add_argument('--min_length', type=int, required=True, help='Minimum transcript length')
    parser.add_argument('--min_exons', type=int, required=True, help='Minimum exon count')
    args = parser.parse_args()

    print("Creating GTF database...")
    db = gffutils.create_db(args.gtf, 'transcripts.db', force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)
    db = gffutils.FeatureDB('transcripts.db')

    filtered_count = 0
    total_count = 0

    with open('filtered_transcripts.gtf', 'w') as out_gtf, open('filtered_transcripts.fasta', 'w') as out_fasta:
        for transcript in db.features_of_type('transcript'):
            total_count += 1
            exons = list(db.children(transcript, featuretype='exon'))
            num_exons = len(exons)
            transcript_length = sum([exon.end - exon.start + 1 for exon in exons])
            if transcript_length >= args.min_length and num_exons >= args.min_exons:
                filtered_count += 1
                out_gtf.write(str(transcript) + '\\n')
                for exon in exons:
                    out_gtf.write(str(exon) + '\\n')
                transcript_id = transcript.attributes.get('transcript_id', [transcript.id])[0]
                seq_record = SeqRecord(Seq('N' * transcript_length), id=transcript_id, description=f"length={transcript_length} exons={num_exons}")
                SeqIO.write(seq_record, out_fasta, 'fasta')

    with open('filter_stats.txt', 'w') as stats:
        stats.write(f"Total transcripts: {total_count}\\n")
        stats.write(f"Filtered transcripts: {filtered_count}\\n")
        stats.write(f"Percentage retained: {100*filtered_count/total_count:.2f}%\\n")
        stats.write(f"Min length: {args.min_length}\\n")
        stats.write(f"Min exons: {args.min_exons}\\n")

    print(f"Filtered {filtered_count} out of {total_count} transcripts")


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.argv = [
            "filter_transcripts.py",
            "--gtf",
            "$gtf",
            "--min_length",
            "$min_length",
            "--min_exons",
            "$min_exons",
        ]
    main()
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
