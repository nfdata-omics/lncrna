#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path
import platform


def parse_featurecounts(file_path):
    """Parse featureCounts output."""
    df = pd.read_csv(file_path, sep='\\t', comment='#')
    # Column names: Geneid, Chr, Start, End, Strand, Length, counts
    sample_name = Path(file_path).stem
    if sample_name.endswith('.featureCounts'):
        sample_name = sample_name.replace('.featureCounts', '')
    counts = df.set_index('Geneid').iloc[:, -1]  # Last column is counts
    counts.name = sample_name
    return counts


def parse_htseq(file_path):
    """Parse HTSeq output."""
    df = pd.read_csv(file_path, sep='\\t', header=None, names=['gene_id', 'count'])
    # Remove HTSeq summary rows
    df = df[~df['gene_id'].str.startswith('__')]
    sample_name = Path(file_path).stem
    counts = df.set_index('gene_id')['count']
    counts.name = sample_name
    return counts


def parse_kallisto(file_path):
    """Parse Kallisto abundance.tsv output."""
    df = pd.read_csv(file_path, sep='\\t')
    # Columns: target_id, length, eff_length, est_counts, tpm
    sample_name = Path(file_path).parent.name
    counts = df.set_index('target_id')['est_counts']
    counts.name = sample_name
    return counts


def parse_salmon(file_path):
    """Parse Salmon quant.sf output."""
    df = pd.read_csv(file_path, sep='\\t')
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

            fields = line.strip().split('\\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            # Allow gene, transcript or exon to ensure we catch gene info
            if feature not in ['gene', 'transcript', 'exon']:
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
                # Prioritize 'gene' feature, or take if not yet present
                if gene_id not in gene_info or feature == 'gene':
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
    counts_str = "$count_files"
    gtf = "$gtf"
    method = "$method"
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "expression"

    count_files = counts_str.split()

    print(f"Merging {len(count_files)} count files using {method} format...", file=sys.stderr)

    count_matrix = merge_counts(count_files, method)

    print("Extracting gene information from GTF...", file=sys.stderr)
    gene_info = parse_gtf_gene_info(gtf)

    output = f"{prefix}.count_matrix.tsv"
    gene_info_path = f"{prefix}.gene_info.tsv"
    summary_path = f"{prefix}.sample_summary.tsv"

    print(f"Writing count matrix to {output}...", file=sys.stderr)
    count_matrix.to_csv(output, sep='\\t')

    print(f"Writing gene info to {gene_info_path}...", file=sys.stderr)
    gene_info.to_csv(gene_info_path, sep='\\t')

    if not gene_info.empty and 'gene_type' in gene_info.columns:
        lncrna_ids = gene_info.index[gene_info['gene_type'] == 'lncRNA']
        if len(lncrna_ids) > 0:
            lncrna_matrix = count_matrix.loc[count_matrix.index.intersection(lncrna_ids)]
            lncrna_output = f"{prefix}.lncrna_matrix.tsv"
            print(f"Writing lncRNA-only count matrix to {lncrna_output}...", file=sys.stderr)
            lncrna_matrix.to_csv(lncrna_output, sep='\\t')

    print(f"Writing summary to {summary_path}...", file=sys.stderr)
    with open(summary_path, 'w') as f:
        f.write("Sample\\tTotal_Counts\\tGenes_Detected\\n")
        for sample in count_matrix.columns:
            total = count_matrix[sample].sum()
            detected = (count_matrix[sample] > 0).sum()
            f.write(f"{sample}\\t{total}\\t{detected}\\n")

    print(f"Done! Count matrix shape: {count_matrix.shape}", file=sys.stderr)


def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


if __name__ == '__main__':
    main()
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
