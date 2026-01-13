#!/usr/bin/env python3
"""
Merge individual count files into a single count matrix for differential expression analysis.
Supports featureCounts, HTSeq, Kallisto, and Salmon formats.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def parse_featurecounts(file_path):
    """Parse featureCounts output."""
    df = pd.read_csv(file_path, sep='\t', comment='#')
    # Column names: Geneid, Chr, Start, End, Strand, Length, counts
    sample_name = Path(file_path).stem
    counts = df.set_index('Geneid').iloc[:, -1]  # Last column is counts
    counts.name = sample_name
    return counts


def parse_htseq(file_path):
    """Parse HTSeq output."""
    df = pd.read_csv(file_path, sep='\t', header=None, names=['gene_id', 'count'])
    # Remove HTSeq summary rows
    df = df[~df['gene_id'].str.startswith('__')]
    sample_name = Path(file_path).stem
    counts = df.set_index('gene_id')['count']
    counts.name = sample_name
    return counts


def parse_kallisto(file_path):
    """Parse Kallisto abundance.tsv output."""
    df = pd.read_csv(file_path, sep='\t')
    # Columns: target_id, length, eff_length, est_counts, tpm
    sample_name = Path(file_path).parent.name
    counts = df.set_index('target_id')['est_counts']
    counts.name = sample_name
    return counts


def parse_salmon(file_path):
    """Parse Salmon quant.sf output."""
    df = pd.read_csv(file_path, sep='\t')
    # Columns: Name, Length, EffectiveLength, TPM, NumReads
    sample_name = Path(file_path).parent.name
    counts = df.set_index('Name')['NumReads']
    counts.name = sample_name
    return counts


def parse_gtf_gene_info(gtf_file):
    """Extract gene information from GTF."""
    gene_info = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature != 'gene':
                continue

            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')

            gene_id = attr_dict.get('gene_id', '')
            gene_name = attr_dict.get('gene_name', gene_id)
            gene_type = attr_dict.get('gene_type') or attr_dict.get('gene_biotype', 'unknown')

            if gene_id:
                gene_info[gene_id] = {
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                }

    return pd.DataFrame.from_dict(gene_info, orient='index')


def merge_counts(count_files, method):
    """Merge multiple count files into a single matrix."""

    count_dfs = []

    for count_file in count_files:
        try:
            if method == 'featurecounts':
                counts = parse_featurecounts(count_file)
            elif method == 'htseq':
                counts = parse_htseq(count_file)
            elif method == 'kallisto':
                counts = parse_kallisto(count_file)
            elif method == 'salmon':
                counts = parse_salmon(count_file)
            else:
                raise ValueError(f"Unknown method: {method}")

            count_dfs.append(counts)
            print(f"Parsed {count_file}: {len(counts)} genes", file=sys.stderr)

        except Exception as e:
            print(f"Warning: Could not parse {count_file}: {e}", file=sys.stderr)
            continue

    if not count_dfs:
        raise ValueError("No valid count files found")

    # Merge all counts into a single dataframe
    count_matrix = pd.concat(count_dfs, axis=1)
    count_matrix = count_matrix.fillna(0).astype(int)

    return count_matrix


def main():
    parser = argparse.ArgumentParser(description='Merge count files into expression matrix')
    parser.add_argument('--counts', nargs='+', required=True, help='Input count files')
    parser.add_argument('--gtf', required=True, help='GTF annotation file')
    parser.add_argument('--method', required=True, choices=['featurecounts', 'htseq', 'kallisto', 'salmon'],
                        help='Quantification method used')
    parser.add_argument('--output', required=True, help='Output count matrix file')
    parser.add_argument('--gene_info', required=True, help='Output gene information file')
    parser.add_argument('--summary', required=True, help='Output summary file')
    args = parser.parse_args()

    print(f"Merging {len(args.counts)} count files using {args.method} format...", file=sys.stderr)

    # Merge count files
    count_matrix = merge_counts(args.counts, args.method)

    # Parse gene information from GTF
    print("Extracting gene information from GTF...", file=sys.stderr)
    gene_info = parse_gtf_gene_info(args.gtf)

    # Write count matrix
    print(f"Writing count matrix to {args.output}...", file=sys.stderr)
    count_matrix.to_csv(args.output, sep='\t')

    # Write gene info
    print(f"Writing gene info to {args.gene_info}...", file=sys.stderr)
    gene_info.to_csv(args.gene_info, sep='\t')

    # Write summary
    print(f"Writing summary to {args.summary}...", file=sys.stderr)
    with open(args.summary, 'w') as f:
        f.write("Sample\tTotal_Counts\tGenes_Detected\n")
        for sample in count_matrix.columns:
            total = count_matrix[sample].sum()
            detected = (count_matrix[sample] > 0).sum()
            f.write(f"{sample}\t{total}\t{detected}\n")

    print(f"Done! Count matrix shape: {count_matrix.shape}", file=sys.stderr)


if __name__ == '__main__':
    main()
