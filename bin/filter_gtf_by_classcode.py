#!/usr/bin/env python3

"""
Filter GTF file by gffcompare class codes
Designed to extract novel lncRNA candidates from gffcompare output
"""

import argparse
import sys
import re
from collections import defaultdict


def parse_classcode_from_attributes(attributes_str):
    """
    Extract class_code from GTF attributes string

    gffcompare adds: class_code "i"; to the attributes
    """
    match = re.search(r'class_code "([^"]+)"', attributes_str)
    if match:
        return match.group(1)
    return None


def extract_attribute(attributes_str, attribute_name):
    """Extract any attribute value from GTF attributes string"""
    pattern = f'{attribute_name} "([^"]+)"'
    match = re.search(pattern, attributes_str)
    if match:
        return match.group(1)
    return None


def filter_gtf_by_classcode(gtf_file, output_file, target_classcodes, stats_file):
    """
    Filter GTF file to keep only transcripts with specified class codes

    Args:
        gtf_file: Input GTF from gffcompare
        output_file: Output filtered GTF
        target_classcodes: Set of class codes to keep (e.g., {'i', 'u', 'x'})
        stats_file: Output statistics file
    """

    # Track statistics
    stats = defaultdict(int)
    stats['total_lines'] = 0
    stats['header_lines'] = 0
    stats['transcript_lines'] = 0
    stats['transcripts_by_classcode'] = defaultdict(int)
    stats['kept_transcripts'] = 0
    stats['filtered_transcripts'] = 0

    # Track which transcript IDs to keep (first pass)
    transcripts_to_keep = set()

    print(f"First pass: Identifying transcripts with class codes: {', '.join(sorted(target_classcodes))}",
          file=sys.stderr)

    # First pass: identify transcripts with target class codes
    with open(gtf_file, 'r') as f:
        for line in f:
            stats['total_lines'] += 1

            # Keep header lines
            if line.startswith('#'):
                stats['header_lines'] += 1
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            # Only process transcript features for classification
            if feature_type == 'transcript':
                stats['transcript_lines'] += 1

                # Extract class code
                class_code = parse_classcode_from_attributes(attributes)

                if class_code:
                    stats['transcripts_by_classcode'][class_code] += 1

                    # Check if this class code should be kept
                    if class_code in target_classcodes:
                        transcript_id = extract_attribute(attributes, 'transcript_id')
                        if transcript_id:
                            transcripts_to_keep.add(transcript_id)
                            stats['kept_transcripts'] += 1
                    else:
                        stats['filtered_transcripts'] += 1
                else:
                    # No class code found - keep it (might be reference annotation)
                    transcript_id = extract_attribute(attributes, 'transcript_id')
                    if transcript_id:
                        transcripts_to_keep.add(transcript_id)

    print(f"Identified {len(transcripts_to_keep)} transcripts to keep", file=sys.stderr)

    # Second pass: write filtered GTF
    print(f"Second pass: Writing filtered GTF to {output_file}", file=sys.stderr)

    written_lines = 0

    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Always keep header lines
            if line.startswith('#'):
                outfile.write(line)
                written_lines += 1
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            # Extract transcript_id
            transcript_id = extract_attribute(attributes, 'transcript_id')

            # Keep line if transcript is in our keep list
            if transcript_id and transcript_id in transcripts_to_keep:
                outfile.write(line)
                written_lines += 1

    print(f"Wrote {written_lines} lines to {output_file}", file=sys.stderr)
    if written_lines == 0:
        print(f"No transcripts matched target class codes; empty output: {output_file}", file=sys.stderr)
        sys.exit(2)

    # Write statistics
    with open(stats_file, 'w') as stats_out:
        stats_out.write("=" * 60 + "\n")
        stats_out.write("GTF FILTERING BY CLASS CODE - STATISTICS\n")
        stats_out.write("=" * 60 + "\n\n")

        stats_out.write(f"Input file: {gtf_file}\n")
        stats_out.write(f"Output file: {output_file}\n")
        stats_out.write(f"Target class codes: {', '.join(sorted(target_classcodes))}\n\n")

        stats_out.write("INPUT STATISTICS:\n")
        stats_out.write("-" * 60 + "\n")
        stats_out.write(f"Total lines: {stats['total_lines']:,}\n")
        stats_out.write(f"Header lines: {stats['header_lines']:,}\n")
        stats_out.write(f"Transcript lines: {stats['transcript_lines']:,}\n\n")

        stats_out.write("TRANSCRIPTS BY CLASS CODE:\n")
        stats_out.write("-" * 60 + "\n")
        for classcode in sorted(stats['transcripts_by_classcode'].keys()):
            count = stats['transcripts_by_classcode'][classcode]
            percentage = (count / stats['transcript_lines'] * 100) if stats['transcript_lines'] > 0 else 0
            kept_marker = " ✓ KEPT" if classcode in target_classcodes else " ✗ FILTERED"
            stats_out.write(f"  {classcode}: {count:,} ({percentage:.2f}%){kept_marker}\n")

        stats_out.write("\n")
        stats_out.write("OUTPUT STATISTICS:\n")
        stats_out.write("-" * 60 + "\n")
        stats_out.write(f"Transcripts kept: {stats['kept_transcripts']:,}\n")
        stats_out.write(f"Transcripts filtered: {stats['filtered_transcripts']:,}\n")
        stats_out.write(f"Lines written: {written_lines:,}\n\n")

        if stats['transcript_lines'] > 0:
            kept_pct = (stats['kept_transcripts'] / stats['transcript_lines'] * 100)
            stats_out.write(f"Retention rate: {kept_pct:.2f}%\n")


def main():
    parser = argparse.ArgumentParser(
        description='Filter GTF file by gffcompare class codes to extract novel lncRNA candidates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Class Code Meanings (gffcompare):
  i = Intronic       : Transcript within intron of reference gene
  u = Intergenic     : Novel intergenic transcript
  x = Antisense      : Transcript on opposite strand of reference gene

  Other codes (not kept by default):
  = = Complete match
  c = Contained in reference
  j = Novel splice junction
  e = Single exon transfrag overlapping reference exon
  o = Other overlap with reference

Examples:
  # Default: Filter for potential lncRNAs (i, u, x)
  filter_gtf_by_classcode.py --gtf merged.annotated.gtf --prefix lncrna_candidates

  # Custom class codes
  filter_gtf_by_classcode.py --gtf merged.annotated.gtf --classcodes i,u --prefix intergenic_only
        """
    )

    parser.add_argument('--gtf', type=str, required=True,
                        help='Input GTF file from gffcompare (merged.annotated.gtf)')
    parser.add_argument('--classcodes', type=str, default='i,u,x',
                        help='Comma-separated class codes to keep (default: i,u,x)')
    parser.add_argument('--prefix', type=str, required=True,
                        help='Output prefix')

    args = parser.parse_args()

    # Parse class codes
    target_classcodes = set(code.strip() for code in args.classcodes.split(','))

    output_gtf = f"{args.prefix}.filtered.gtf"
    stats_file = f"{args.prefix}.classcode_stats.txt"

    try:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"Filtering GTF by class codes: {', '.join(sorted(target_classcodes))}", file=sys.stderr)
        print(f"{'='*60}\n", file=sys.stderr)

        filter_gtf_by_classcode(args.gtf, output_gtf, target_classcodes, stats_file)

        print(f"\n✓ Successfully created filtered GTF: {output_gtf}", file=sys.stderr)
        print(f"✓ Statistics saved to: {stats_file}", file=sys.stderr)

    except Exception as e:
        print(f"\n✗ Error filtering GTF: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
