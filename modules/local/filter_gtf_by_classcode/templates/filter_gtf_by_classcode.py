#!/usr/bin/env python3

"""
Filter GTF file by gffcompare class codes
Designed to extract novel lncRNA candidates from gffcompare output

For class code 'j' (novel junction): additionally filters out transcripts
whose reference gene (ref_gene_id) is protein_coding in the reference GTF.
"""

import argparse
import sys
import re
from collections import defaultdict
import platform


def parse_classcode_from_attributes(attributes_str):
    match = re.search(r'class_code "([^"]+)"', attributes_str)
    if match:
        return match.group(1)
    return None


def extract_attribute(attributes_str, attribute_name):
    pattern = f'{attribute_name} "([^"]+)"'
    match = re.search(pattern, attributes_str)
    if match:
        return match.group(1)
    return None


def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def load_protein_coding_genes(ref_gtf):
    """
    Load set of protein_coding gene IDs from reference GTF.
    Used to filter out 'j' class code transcripts overlapping protein_coding genes.
    """
    protein_coding = set()
    with open(ref_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue
            attr = fields[8]
            biotype = (
                re.search(r'gene_biotype "([^"]+)"', attr) or
                re.search(r'gene_type "([^"]+)"', attr)
            )
            if biotype and biotype.group(1) == 'protein_coding':
                gene_id = re.search(r'gene_id "([^"]+)"', attr)
                if gene_id:
                    protein_coding.add(gene_id.group(1))

    print(f"Loaded {len(protein_coding)} protein_coding genes from reference GTF", file=sys.stderr)
    return protein_coding


def filter_gtf_by_classcode(gtf_file, output_file, target_classcodes, stats_file, protein_coding_genes=None):
    """
    Filter GTF file to keep only transcripts with specified class codes.

    For class code 'j': additionally excludes transcripts whose ref_gene_id
    is protein_coding in the reference GTF (requires protein_coding_genes set).
    """

    stats = defaultdict(int)
    stats['total_lines'] = 0
    stats['header_lines'] = 0
    stats['transcript_lines'] = 0
    stats['transcripts_by_classcode'] = defaultdict(int)
    stats['kept_transcripts'] = 0
    stats['filtered_transcripts'] = 0
    stats['j_filtered_protein_coding'] = 0

    transcripts_to_keep = set()

    print(f"First pass: Identifying transcripts with class codes: {', '.join(sorted(target_classcodes))}",
          file=sys.stderr)

    with open(gtf_file, 'r') as f:
        for line in f:
            stats['total_lines'] += 1

            if line.startswith('#'):
                stats['header_lines'] += 1
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            if feature_type == 'transcript':
                stats['transcript_lines'] += 1
                class_code = parse_classcode_from_attributes(attributes)

                if class_code:
                    stats['transcripts_by_classcode'][class_code] += 1

                    if class_code not in target_classcodes:
                        stats['filtered_transcripts'] += 1
                        continue

                    # For class code 'j': filter out protein_coding reference genes
                    if class_code == 'j' and protein_coding_genes is not None:
                        ref_gene_id = extract_attribute(attributes, 'ref_gene_id')
                        if ref_gene_id and ref_gene_id in protein_coding_genes:
                            stats['j_filtered_protein_coding'] += 1
                            stats['filtered_transcripts'] += 1
                            continue

                    transcript_id = extract_attribute(attributes, 'transcript_id')
                    if transcript_id:
                        transcripts_to_keep.add(transcript_id)
                        stats['kept_transcripts'] += 1

                else:
                    # No class code — keep (reference annotation lines)
                    transcript_id = extract_attribute(attributes, 'transcript_id')
                    if transcript_id:
                        transcripts_to_keep.add(transcript_id)

    print(f"Identified {len(transcripts_to_keep)} transcripts to keep", file=sys.stderr)

    # Second pass: write filtered GTF
    print(f"Second pass: Writing filtered GTF to {output_file}", file=sys.stderr)
    written_lines = 0

    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                written_lines += 1
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            attributes = fields[8]
            transcript_id = extract_attribute(attributes, 'transcript_id')

            if transcript_id and transcript_id in transcripts_to_keep:
                outfile.write(line)
                written_lines += 1

    print(f"Wrote {written_lines} lines to {output_file}", file=sys.stderr)
    if written_lines == 0:
        print(f"No transcripts matched target class codes; empty output: {output_file}", file=sys.stderr)
        sys.exit(2)

    # Write statistics
    with open(stats_file, 'w') as stats_out:
        stats_out.write("=" * 60 + "\\n")
        stats_out.write("GTF FILTERING BY CLASS CODE - STATISTICS\\n")
        stats_out.write("=" * 60 + "\\n\\n")

        stats_out.write(f"Input file: {gtf_file}\\n")
        stats_out.write(f"Output file: {output_file}\\n")
        stats_out.write(f"Target class codes: {', '.join(sorted(target_classcodes))}\\n\\n")

        stats_out.write("INPUT STATISTICS:\\n")
        stats_out.write("-" * 60 + "\\n")
        stats_out.write(f"Total lines: {stats['total_lines']:,}\\n")
        stats_out.write(f"Header lines: {stats['header_lines']:,}\\n")
        stats_out.write(f"Transcript lines: {stats['transcript_lines']:,}\\n\\n")

        stats_out.write("TRANSCRIPTS BY CLASS CODE:\\n")
        stats_out.write("-" * 60 + "\\n")
        for classcode in sorted(stats['transcripts_by_classcode'].keys()):
            count = stats['transcripts_by_classcode'][classcode]
            percentage = (count / stats['transcript_lines'] * 100) if stats['transcript_lines'] > 0 else 0
            kept_marker = "KEPT" if classcode in target_classcodes else "FILTERED"
            stats_out.write(f"  {classcode}: {count:,} ({percentage:.2f}%){kept_marker}\\n")

        stats_out.write("\\n")
        stats_out.write("OUTPUT STATISTICS:\\n")
        stats_out.write("-" * 60 + "\\n")
        stats_out.write(f"Transcripts kept: {stats['kept_transcripts']:,}\\n")
        stats_out.write(f"Transcripts filtered: {stats['filtered_transcripts']:,}\\n")
        if 'j' in target_classcodes:
            stats_out.write(f"  of which 'j' filtered (protein_coding ref): {stats['j_filtered_protein_coding']:,}\\n")
        stats_out.write(f"Lines written: {written_lines:,}\\n\\n")

        if stats['transcript_lines'] > 0:
            kept_pct = (stats['kept_transcripts'] / stats['transcript_lines'] * 100)
            stats_out.write(f"Retention rate: {kept_pct:.2f}%\\n")


def main():
    parser = argparse.ArgumentParser(
        description='Filter GTF file by gffcompare class codes to extract novel lncRNA candidates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Class Code Meanings (gffcompare):
  i = Intronic       : Fully contained in a reference intron
  u = Intergenic     : Novel intergenic transcript
  x = Antisense      : Exonic overlap on opposite strand
  j = Novel junction : At least one junction match with reference gene
                       (filtered for protein_coding if --ref_gtf provided)
        """
    )

    parser.add_argument('--gtf',     type=str, required=True,
                        help='Input GTF file from gffcompare (merged.annotated.gtf)')
    parser.add_argument('--prefix',  type=str, required=True,
                        help='Output prefix')
    parser.add_argument('--ref_gtf', type=str, required=False, default=None,
                        help='Reference annotation GTF. Used to filter j class code '
                             'transcripts overlapping protein_coding genes.')

    args = parser.parse_args()

    classcodes = 'i,u,x,j'
    target_classcodes = set(code.strip() for code in classcodes.split(','))

    output_gtf = f"{args.prefix}.filtered.gtf"
    stats_file = f"{args.prefix}.classcode_stats.txt"

    # Load protein_coding genes if ref_gtf provided
    protein_coding_genes = None
    if args.ref_gtf:
        protein_coding_genes = load_protein_coding_genes(args.ref_gtf)

    try:
        print(f"\\n{'='*60}", file=sys.stderr)
        print(f"Filtering GTF by class codes: {', '.join(sorted(target_classcodes))}", file=sys.stderr)
        if protein_coding_genes:
            print(f"Filtering 'j' transcripts overlapping protein_coding genes", file=sys.stderr)
        print(f"{'='*60}\\n", file=sys.stderr)

        filter_gtf_by_classcode(args.gtf, output_gtf, target_classcodes, stats_file, protein_coding_genes)

        print(f"Successfully created filtered GTF: {output_gtf}", file=sys.stderr)
        print(f"Statistics saved to: {stats_file}", file=sys.stderr)

    except Exception as e:
        print(f"Error filtering GTF: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.argv = [
            "filter_gtf_by_classcode.py",
            "--gtf",    "$gtf",
            "--prefix", "${task.ext.prefix}" if "${task.ext.prefix}" != "null" else "${gtf.baseName}",
            "--ref_gtf","$ref_gtf",
        ]
    main()
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
