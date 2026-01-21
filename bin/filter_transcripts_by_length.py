#!/usr/bin/env python3

"""
Filter transcripts by minimum length (default: 200 nt)
lncRNAs are typically >= 200 nucleotides
"""

import argparse
import sys
from Bio import SeqIO


def filter_transcripts_by_length(gtf_file, fasta_file, min_length, output_prefix):
    """Filter transcripts shorter than min_length"""

    # Read FASTA and identify transcripts that pass length filter
    passing_transcripts = set()
    stats = {'total': 0, 'passed': 0, 'filtered': 0}

    print(f"Reading sequences from {fasta_file}...", file=sys.stderr)

    with open(f"{output_prefix}.length_filtered.fa", 'w') as fasta_out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            stats['total'] += 1
            seq_length = len(record.seq)

            if seq_length >= min_length:
                SeqIO.write(record, fasta_out, "fasta")
                passing_transcripts.add(record.id)
                stats['passed'] += 1
            else:
                stats['filtered'] += 1

    print(f"Filtered {stats['passed']} / {stats['total']} transcripts (>= {min_length} nt)",
          file=sys.stderr)

    if stats['passed'] == 0:
        with open(f"{output_prefix}.length_filtered.gtf", 'w') as gtf_out:
            pass
        with open(f"{output_prefix}.length_filtered.fa", 'w') as fasta_out:
            pass
        with open(f"{output_prefix}.length_stats.txt", 'w') as stats_out:
            stats_out.write("LENGTH FILTERING STATISTICS\n")
            stats_out.write("=" * 50 + "\n\n")
            stats_out.write(f"Minimum length threshold: {min_length} nt\n\n")
            stats_out.write(f"Total transcripts: {stats['total']}\n")
            stats_out.write(f"Passed filter: 0 (0.00%)\n")
            stats_out.write(f"Filtered out: {stats['filtered']} (100.00%)\n")
        print("ERROR: No transcripts pass the length filter", file=sys.stderr)
        sys.exit(2)

    # Filter GTF based on passing transcripts
    print(f"Filtering GTF file...", file=sys.stderr)

    with open(gtf_file, 'r') as gtf_in, \
         open(f"{output_prefix}.length_filtered.gtf", 'w') as gtf_out:

        for line in gtf_in:
            if line.startswith('#'):
                gtf_out.write(line)
                continue

            # Extract transcript_id
            if 'transcript_id' in line:
                import re
                match = re.search(r'transcript_id "([^"]+)"', line)
                if match:
                    transcript_id = match.group(1)
                    if transcript_id in passing_transcripts:
                        gtf_out.write(line)

    # Write statistics
    with open(f"{output_prefix}.length_stats.txt", 'w') as stats_out:
        stats_out.write("LENGTH FILTERING STATISTICS\n")
        stats_out.write("=" * 50 + "\n\n")
        stats_out.write(f"Minimum length threshold: {min_length} nt\n\n")
        stats_out.write(f"Total transcripts: {stats['total']}\n")
        pct_passed = (stats['passed'] / stats['total'] * 100) if stats['total'] > 0 else 0
        pct_filtered = (stats['filtered'] / stats['total'] * 100) if stats['total'] > 0 else 0
        stats_out.write(f"Passed filter: {stats['passed']} ({pct_passed:.2f}%)\n")
        stats_out.write(f"Filtered out: {stats['filtered']} ({pct_filtered:.2f}%)\n")

    print(f"✓ Done: {output_prefix}.length_filtered.gtf", file=sys.stderr)
    print(f"✓ Done: {output_prefix}.length_filtered.fa", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Filter transcripts by minimum length (lncRNAs >= 200 nt)'
    )

    parser.add_argument('--gtf', required=True, help='Input GTF file')
    parser.add_argument('--fasta', required=True, help='Input transcript FASTA')
    parser.add_argument('--min_length', type=int, default=200,
                        help='Minimum transcript length (default: 200)')
    parser.add_argument('--prefix', required=True, help='Output prefix')

    args = parser.parse_args()

    filter_transcripts_by_length(args.gtf, args.fasta, args.min_length, args.prefix)


if __name__ == '__main__':
    main()
