#!/usr/bin/env python3


"""
Analyze Trans relationships between lncRNAs and protein-coding genes.
Calculates expression correlation (Pearson/Spearman) for genes on different
chromosomes or beyond a defined cis window.
"""


import argparse
import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
import platform
try:
    from statsmodels.stats.multitest import multipletests
    STATSMODELS_AVAIL = True
except ImportError:
    STATSMODELS_AVAIL = False


def parse_gtf_coords(gtf_file):
    """
    Parse GTF to get gene coordinates and types.
    Returns: dict {gene_id: {'chrom': str, 'start': int, 'end': int, 'type': str, 'name': str}}
    """
    genes = {}
    print(f"Parsing GTF: {gtf_file}...", file=sys.stderr)

    lines_read = 0
    features_found = set()

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            lines_read += 1
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            features_found.add(feature)

            if feature not in ['gene', 'transcript']:
                continue

            chrom = fields[0]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue

            attributes = fields[8]

            attr_dict = {}
            for attr in attributes.strip().rstrip(';').split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')

            gene_id = attr_dict.get('gene_id')
            if not gene_id:
                continue

            gene_name = attr_dict.get('gene_name', gene_id)
            gene_type = (attr_dict.get('gene_type') or
                         attr_dict.get('gene_biotype') or
                         attr_dict.get('biotype') or
                         attr_dict.get('transcript_type') or
                         attr_dict.get('transcript_biotype'))

            if feature == 'gene':
                genes[gene_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'type': gene_type,
                    'name': gene_name
                }
            elif feature == 'transcript':
                if gene_id not in genes:
                    genes[gene_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'type': gene_type,
                        'name': gene_name
                    }
                else:
                    genes[gene_id]['start'] = min(genes[gene_id]['start'], start)
                    genes[gene_id]['end'] = max(genes[gene_id]['end'], end)
                    if not genes[gene_id].get('type') and gene_type:
                        genes[gene_id]['type'] = gene_type

    print(f"DEBUG: Parsed {lines_read} lines. Features found: {features_found}", file=sys.stderr)
    return genes


def calculate_correlation(s1, s2, method='pearson'):
    if method == 'pearson':
        return stats.pearsonr(s1, s2)
    elif method == 'spearman':
        return stats.spearmanr(s1, s2)
    else:
        raise ValueError(f"Unknown method: {method}")


def is_trans(lnc_info, coding_info, cis_window):
    """
    A pair is Trans if they are on different chromosomes,
    or on the same chromosome but beyond the cis window.
    """
    if lnc_info['chrom'] != coding_info['chrom']:
        return True

    start1, end1 = lnc_info['start'], lnc_info['end']
    start2, end2 = coding_info['start'], coding_info['end']

    if end1 < start2:
        dist = start2 - end1
    elif end2 < start1:
        dist = start1 - end2
    else:
        dist = 0  # Overlap

    return dist > cis_window


def main():
    expression = "$expression_matrix"
    gtf = "$gtf"
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "${meta.id}"
    output_trans = f"{prefix}.trans_results.tsv"

    extra_args = "$task.ext.args"

    if len(sys.argv) == 1:
        import shlex
        argv = [
            sys.argv[0],
            "--expression", expression,
            "--gtf", gtf,
            "--output_trans", output_trans,
        ]
        if extra_args and extra_args != "null":
            argv += shlex.split(extra_args)
        sys.argv = argv

    parser = argparse.ArgumentParser(description='Analyze Trans lncRNA-gene correlations')
    parser.add_argument('--expression', required=True, help='Expression matrix (genes x samples)')
    parser.add_argument('--gtf', required=True, help='GTF file for coordinates')
    parser.add_argument('--output_trans', required=True, help='Output file for Trans results')
    parser.add_argument('--cis_window', type=int, default=100000, help='Cis window size to exclude (bp) [default: 100000]')
    parser.add_argument('--method', choices=['pearson', 'spearman'], default='pearson', help='Correlation method')
    parser.add_argument('--min_corr', type=float, default=0.3, help='Minimum absolute correlation to report')
    parser.add_argument('--pvalue', type=float, default=0.05, help='P-value threshold')
    parser.add_argument('--chunk_size', type=int, default=50, help='Number of lncRNAs to process per chunk [default: 50]')

    args = parser.parse_args()

    # 1. Load Expression Data
    print(f"Loading expression matrix: {args.expression}...", file=sys.stderr)
    df = pd.read_csv(args.expression, sep='\t', index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').dropna()

    if df.shape[1] > 1:
        variance = df.var(axis=1)
        df = df.loc[variance > 0]
        print(f"DEBUG: Retained {len(df)} genes with non-zero variance.", file=sys.stderr)

    if df.shape[1] < 2:
        print("Warning: Expression matrix has fewer than 2 samples. Correlation analysis cannot be performed.", file=sys.stderr)
        pd.DataFrame(columns=['lncRNA_ID', 'lncRNA_Name', 'Gene_ID', 'Gene_Name',
                               'lncRNA_Chrom', 'Gene_Chrom', 'Correlation', 'P_value', 'FDR']).to_csv(args.output_trans, sep='\t', index=False)
        return

    if df.shape[1] < 3:
        print(f"Warning: Expression matrix has only {df.shape[1]} samples. P-values will not be statistically meaningful.", file=sys.stderr)

    # 2. Load Gene Annotations
    gene_info = parse_gtf_coords(args.gtf)

    print(f"DEBUG: First 5 expression IDs: {list(df.index[:5])}", file=sys.stderr)
    print(f"DEBUG: First 5 GTF IDs: {list(gene_info.keys())[:5]}", file=sys.stderr)

    # 3. Separate lncRNA and Protein Coding
    lncrna_ids = []
    coding_ids = []

    def clean_id(x):
        return str(x).split('.')[0]

    expr_id_map = {clean_id(gid): gid for gid in df.index}
    valid_genes_clean = set(expr_id_map.keys())

    for gid, info in gene_info.items():
        gid_clean = clean_id(gid)

        if gid_clean not in valid_genes_clean:
            continue

        original_expr_id = expr_id_map[gid_clean]
        gtype = info.get('type', '')

        if gtype in ['lncRNA', 'novel_lncRNA', 'lincRNA', 'antisense', 'processed_transcript',
                     'sense_intronic', 'sense_overlapping', '3prime_overlapping_ncRNA']:
            lncrna_ids.append(original_expr_id)
        elif gtype in ['protein_coding', 'mRNA']:
            coding_ids.append(original_expr_id)

    print(f"Found {len(lncrna_ids)} lncRNAs and {len(coding_ids)} protein-coding genes in expression matrix.", file=sys.stderr)

    if not lncrna_ids or not coding_ids:
        print("Error: Not enough genes found to analyze.", file=sys.stderr)
        types_found = {}
        for gid, info in gene_info.items():
            t = info.get('type', 'unknown')
            types_found[t] = types_found.get(t, 0) + 1
        print(f"DEBUG: Gene types found in GTF: {types_found}", file=sys.stderr)
        sys.exit(1)

    # 4. Trans Analysis
    print(f"Running Trans analysis (Cis exclusion window: {args.cis_window} bp)...", file=sys.stderr)
    print(f"Total pairs to evaluate: {len(lncrna_ids)} lncRNAs x {len(coding_ids)} coding genes = {len(lncrna_ids) * len(coding_ids):,}", file=sys.stderr)

    trans_results = []
    total_lnc = len(lncrna_ids)

    # Process in chunks to manage memory
    for chunk_start in range(0, total_lnc, args.chunk_size):
        chunk = lncrna_ids[chunk_start:chunk_start + args.chunk_size]
        print(f"Processing lncRNAs {chunk_start + 1}-{min(chunk_start + args.chunk_size, total_lnc)} of {total_lnc}...", file=sys.stderr)

        for lnc in chunk:
            lnc_info = gene_info[lnc]

            for cid in coding_ids:
                coding_info = gene_info[cid]

                # Only keep trans pairs
                if not is_trans(lnc_info, coding_info, args.cis_window):
                    continue

                r, p = calculate_correlation(df.loc[lnc], df.loc[cid], args.method)

                if np.isnan(r):
                    continue

                is_significant = (p <= args.pvalue)
                if df.shape[1] < 3:
                    is_significant = True

                if abs(r) >= args.min_corr and is_significant:
                    trans_results.append({
                        'lncRNA_ID': lnc,
                        'lncRNA_Name': lnc_info['name'],
                        'Gene_ID': cid,
                        'Gene_Name': coding_info['name'],
                        'lncRNA_Chrom': lnc_info['chrom'],
                        'Gene_Chrom': coding_info['chrom'],
                        'Correlation': round(r, 6),
                        'P_value': p
                    })

    # 5. Build DataFrame and apply FDR
    trans_df = pd.DataFrame(trans_results)

    if not trans_df.empty:
        if STATSMODELS_AVAIL and df.shape[1] >= 3:
            try:
                reject, qvals, _, _ = multipletests(trans_df['P_value'], method='fdr_bh')
                trans_df['FDR'] = qvals
            except Exception as e:
                print(f"Warning: FDR correction failed: {e}", file=sys.stderr)
                trans_df['FDR'] = trans_df['P_value']
        else:
            trans_df['FDR'] = trans_df['P_value']
            print("Warning: statsmodels not available, skipping FDR correction.", file=sys.stderr)

        trans_df = trans_df.sort_values('FDR')
    else:
        trans_df = pd.DataFrame(columns=['lncRNA_ID', 'lncRNA_Name', 'Gene_ID', 'Gene_Name',
                                         'lncRNA_Chrom', 'Gene_Chrom', 'Correlation', 'P_value', 'FDR'])

    print(f"Writing {len(trans_df)} trans-interactions to {args.output_trans}", file=sys.stderr)
    trans_df.to_csv(args.output_trans, sep='\t', index=False)


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
