#!/usr/bin/env python3
"""
Filter lncRNAs based on BLAST results against protein databases.
Remove transcripts with significant similarity to protein-coding genes.
"""

import argparse
from Bio import SeqIO
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


def parse_blast_hits(blast_file):
    """Parse BLAST output and return set of query IDs with hits."""
    hits = set()

    if not blast_file:
        return hits

    try:
        with open(blast_file, 'r') as f:
            for line in f:
                if line.strip():
                    qseqid = line.split('\t')[0]
                    hits.add(qseqid)
    except FileNotFoundError:
        pass

    return hits


def filter_sequences(input_fasta, blast_hits, output_fasta, stats_file):
    """Filter sequences removing those with BLAST hits."""

    total = 0
    kept = 0
    removed = 0

    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            total += 1

            if record.id not in blast_hits:
                SeqIO.write(record, out_f, 'fasta')
                kept += 1
            else:
                removed += 1

    # Write statistics
    with open(stats_file, 'w') as f:
        f.write("=== BLAST Filtering Statistics ===\\n")
        f.write(f"Total input sequences: {total}\\n")
        f.write(f"Sequences with protein homology (removed): {removed}\\n")
        f.write(f"Sequences without hits (kept): {kept}\\n")
        f.write(f"Retention rate: {kept/total*100:.2f}%\\n")

    print(f"Filtered {total} sequences:", file=sys.stderr)
    print(f"  Kept: {kept} ({kept/total*100:.1f}%)", file=sys.stderr)
    print(f"  Removed: {removed} ({removed/total*100:.1f}%)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Filter lncRNAs based on BLAST results')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--blast_hits', required=True, help='BLAST results file')
    parser.add_argument('--output', required=True, help='Output filtered FASTA')
    parser.add_argument('--stats', required=True, help='Output statistics file')
    args = parser.parse_args()

    print("Parsing BLAST hits...", file=sys.stderr)
    blast_hits = parse_blast_hits(args.blast_hits)
    print(f"Found {len(blast_hits)} sequences with protein homology", file=sys.stderr)

    print("Filtering sequences...", file=sys.stderr)
    filter_sequences(args.input, blast_hits, args.output, args.stats)

    print("Done!", file=sys.stderr)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        if "${task.ext.prefix}" != "null":
            prefix_value = "${task.ext.prefix}"
        elif "$meta.id" != "null":
            prefix_value = "$meta.id"
        else:
            prefix_value = "blast_filter"
        sys.argv = [
            "filter_blast_results.py",
            "--input",
            "$lncrna_fasta",
            "--blast_hits",
            "$blast_hits",
            "--output",
            f"{prefix_value}.filtered.fa",
            "--stats",
            f"{prefix_value}.blast_filter_stats.txt",
        ]
    main()
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
