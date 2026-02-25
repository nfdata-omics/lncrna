#!/usr/bin/env python3

"""
Filter single-exon transcripts
Multi-exon lncRNAs are more reliable than single-exon
"""

import argparse
import sys
from collections import defaultdict
from Bio import SeqIO
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


def filter_transcripts_by_exons(gtf_file, fasta_file, min_exons, output_prefix):
    """Filter transcripts with fewer than min_exons"""

    # Count exons per transcript
    exon_counts = defaultdict(int)

    print(f"Counting exons per transcript from {gtf_file}...", file=sys.stderr)

    with open(gtf_file, 'r') as gtf_in:
        for line in gtf_in:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type == 'exon':
                # Extract transcript_id
                import re
                match = re.search(r'transcript_id "([^"]+)"', fields[8])
                if match:
                    transcript_id = match.group(1)
                    exon_counts[transcript_id] += 1

    # Identify passing transcripts
    passing_transcripts = {
        tid for tid, count in exon_counts.items()
        if count >= min_exons
    }

    stats = {
        'total': len(exon_counts),
        'passed': len(passing_transcripts),
        'filtered': len(exon_counts) - len(passing_transcripts)
    }

    print(f"Filtered {stats['passed']} / {stats['total']} transcripts (>= {min_exons} exons)",
          file=sys.stderr)

    if stats['passed'] == 0:
        with open(f"{output_prefix}.exon_filtered.gtf", 'w') as gtf_out:
            pass
        with open(f"{output_prefix}.exon_filtered.fa", 'w') as fasta_out:
            pass
        with open(f"{output_prefix}.exon_stats.txt", 'w') as stats_out:
            stats_out.write("EXON FILTERING STATISTICS\\n")
            stats_out.write("=" * 50 + "\\n\\n")
            stats_out.write(f"Minimum exon threshold: {min_exons}\\n\\n")
            stats_out.write(f"Total transcripts: {stats['total']}\\n")
            stats_out.write("Passed filter: 0 (0.00%)\\n")
            stats_out.write(f"Filtered out (single-exon): {stats['filtered']} (100.00%)\\n")
        print("ERROR: No transcripts pass the exon filter", file=sys.stderr)
        sys.exit(2)

    # Filter FASTA
    print(f"Filtering FASTA file...", file=sys.stderr)

    with open(f"{output_prefix}.exon_filtered.fa", 'w') as fasta_out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in passing_transcripts:
                SeqIO.write(record, fasta_out, "fasta")

    # Filter GTF
    print(f"Filtering GTF file...", file=sys.stderr)

    with open(gtf_file, 'r') as gtf_in, \
         open(f"{output_prefix}.exon_filtered.gtf", 'w') as gtf_out:

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
    with open(f"{output_prefix}.exon_stats.txt", 'w') as stats_out:
        stats_out.write("EXON FILTERING STATISTICS\\n")
        stats_out.write("=" * 50 + "\\n\\n")
        stats_out.write(f"Minimum exon threshold: {min_exons}\\n\\n")
        stats_out.write(f"Total transcripts: {stats['total']}\\n")
        pct_passed = (stats['passed'] / stats['total'] * 100) if stats['total'] > 0 else 0
        pct_filtered = (stats['filtered'] / stats['total'] * 100) if stats['total'] > 0 else 0
        stats_out.write(f"Passed filter: {stats['passed']} ({pct_passed:.2f}%)\\n")
        stats_out.write(f"Filtered out (single-exon): {stats['filtered']} ({pct_filtered:.2f}%)\\n")

    print(f"✓ Done: {output_prefix}.exon_filtered.gtf", file=sys.stderr)
    print(f"✓ Done: {output_prefix}.exon_filtered.fa", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Filter single-exon transcripts (keep multi-exon only)'
    )

    parser.add_argument('--gtf', required=True, help='Input GTF file')
    parser.add_argument('--fasta', required=True, help='Input transcript FASTA')
    parser.add_argument('--min_exons', type=int, default=2,
                        help='Minimum number of exons (default: 2)')
    parser.add_argument('--prefix', required=True, help='Output prefix')

    args = parser.parse_args()

    filter_transcripts_by_exons(args.gtf, args.fasta, args.min_exons, args.prefix)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.argv = [
            "filter_transcripts_by_exons.py",
            "--gtf",
            "$gtf",
            "--fasta",
            "$fasta",
            "--min_exons",
            "$task.ext.args" if "$task.ext.args" != "null" else "2",
            "--prefix",
            "${task.ext.prefix}" if "${task.ext.prefix}" != "null" else "${meta.id}",
        ]
    main()
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
