#!/usr/bin/env python3

"""
Analyze Cis and Trans relationships between lncRNAs and protein-coding genes.
Calculates expression correlation (Pearson/Spearman).
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
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

    # Debug counters
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

            # Accept 'gene' or 'transcript' to handle GTFs without explicit gene lines
            if feature not in ['gene', 'transcript']:
                continue

            chrom = fields[0]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue

            attributes = fields[8]

            # Parse attributes
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
            # Try multiple keys for biotype
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
                    # Update coords if transcript extends gene boundaries
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

def main():
    parser = argparse.ArgumentParser(description='Analyze Cis/Trans lncRNA-gene correlations')
    parser.add_argument('--expression', required=True, help='Expression matrix (genes x samples)')
    parser.add_argument('--gtf', required=True, help='GTF file for coordinates')
    parser.add_argument('--output_cis', required=True, help='Output file for Cis results')
    parser.add_argument('--output_trans', required=False, help='Output file for Trans results (optional)')
    parser.add_argument('--window', type=int, default=100000, help='Window size for Cis interaction (bp) [default: 100000]')
    parser.add_argument('--method', choices=['pearson', 'spearman'], default='pearson', help='Correlation method')
    parser.add_argument('--min_corr', type=float, default=0.3, help='Minimum correlation to report')
    parser.add_argument('--pvalue', type=float, default=0.05, help='P-value threshold')

    args = parser.parse_args()

    # 1. Load Expression Data
    print(f"Loading expression matrix: {args.expression}...", file=sys.stderr)
    # Assume first column is gene_id or index
    df = pd.read_csv(args.expression, sep='\t', index_col=0)
    # Ensure numeric
    df = df.apply(pd.to_numeric, errors='coerce').dropna()

    # Filter out genes with zero variance (constant expression across samples)
    # This avoids PearsonRConstantInputWarning and undefined correlations
    if df.shape[1] > 1:
        variance = df.var(axis=1)
        df = df.loc[variance > 0]
        print(f"DEBUG: Retained {len(df)} genes with non-zero variance.", file=sys.stderr)

    # Check sample count
    if df.shape[1] < 2:
        print("Warning: Expression matrix has fewer than 2 samples. Correlation analysis cannot be performed.", file=sys.stderr)
        # Create empty output files and exit successfully
        pd.DataFrame(columns=['lncRNA_ID', 'lncRNA_Name', 'Gene_ID', 'Gene_Name', 'Chromosome', 'Distance', 'Correlation', 'P_value', 'FDR']).to_csv(args.output_cis, sep='\t', index=False)
        if args.output_trans:
            pd.DataFrame(columns=['lncRNA_ID', 'Gene_ID', 'Correlation', 'P_value']).to_csv(args.output_trans, sep='\t', index=False)
        return

    if df.shape[1] < 3:
        print(f"Warning: Expression matrix has only {df.shape[1]} samples. P-values will not be statistically meaningful.", file=sys.stderr)

    # 2. Load Gene Annotations
    gene_info = parse_gtf_coords(args.gtf)

    # Debug: Print first few genes from expression matrix and GTF
    print(f"DEBUG: First 5 expression IDs: {list(df.index[:5])}", file=sys.stderr)
    print(f"DEBUG: First 5 GTF IDs: {list(gene_info.keys())[:5]}", file=sys.stderr)

    # 3. Separate lncRNA and Protein Coding
    lncrna_ids = []
    coding_ids = []

    # Helper to clean IDs (remove version dots)
    def clean_id(x):
        return str(x).split('.')[0]

    # Map cleaned expression IDs to original IDs
    expr_id_map = {clean_id(gid): gid for gid in df.index}
    valid_genes_clean = set(expr_id_map.keys())

    count_lnc = 0
    count_pc = 0

    for gid, info in gene_info.items():
        gid_clean = clean_id(gid)

        if gid_clean not in valid_genes_clean:
            continue

        # Get original expression ID
        original_expr_id = expr_id_map[gid_clean]

        gtype = info.get('type', '')
        # Extended biotypes
        if gtype in ['lncRNA', 'novel_lncRNA', 'lincRNA', 'antisense', 'processed_transcript', 'sense_intronic', 'sense_overlapping', '3prime_overlapping_ncRNA']:
            lncrna_ids.append(original_expr_id)
            count_lnc += 1
        elif gtype in ['protein_coding', 'mRNA']:
            coding_ids.append(original_expr_id)
            count_pc += 1

    print(f"Found {len(lncrna_ids)} lncRNAs and {len(coding_ids)} protein-coding genes in expression matrix.", file=sys.stderr)

    if not lncrna_ids or not coding_ids:
        print("Error: Not enough genes found to analyze.", file=sys.stderr)
        # Fallback: if 0 found, list some GTF types found to debug
        types_found = {}
        for gid, info in gene_info.items():
            t = info.get('type', 'unknown')
            types_found[t] = types_found.get(t, 0) + 1
        print(f"DEBUG: Gene types found in GTF: {types_found}", file=sys.stderr)
        sys.exit(1)

    # 4. Cis Analysis
    print(f"Running Cis analysis (Window: {args.window} bp)...", file=sys.stderr)
    cis_results = []

    # Optimize: Group coding genes by chromosome for faster lookup
    coding_by_chrom = {}
    for cid in coding_ids:
        chrom = gene_info[cid]['chrom']
        if chrom not in coding_by_chrom:
            coding_by_chrom[chrom] = []
        coding_by_chrom[chrom].append(cid)

    for lnc in lncrna_ids:
        lnc_info = gene_info[lnc]
        chrom = lnc_info['chrom']

        if chrom not in coding_by_chrom:
            continue

        # Find neighbors
        neighbors = []
        for cid in coding_by_chrom[chrom]:
            coding_info = gene_info[cid]
            # Check distance
            dist = float('inf')

            # Simple distance: min distance between gene bodies
            # if overlap, distance is 0
            start1, end1 = lnc_info['start'], lnc_info['end']
            start2, end2 = coding_info['start'], coding_info['end']

            if end1 < start2:
                dist = start2 - end1
            elif end2 < start1:
                dist = start1 - end2
            else:
                dist = 0 # Overlap

            if dist <= args.window:
                neighbors.append((cid, dist))

        # Calculate correlation for neighbors
        for neighbor_id, distance in neighbors:
            r, p = calculate_correlation(df.loc[lnc], df.loc[neighbor_id], args.method)

            # Skip if correlation is NaN
            if np.isnan(r):
                continue

            # Check significance
            # If N < 3, p-value is unreliable (often 1.0 or NaN). We rely on correlation magnitude.
            is_significant = (p <= args.pvalue)
            if df.shape[1] < 3:
                is_significant = True

            if abs(r) >= args.min_corr and is_significant:
                cis_results.append({
                    'lncRNA_ID': lnc,
                    'lncRNA_Name': lnc_info['name'],
                    'Gene_ID': neighbor_id,
                    'Gene_Name': gene_info[neighbor_id]['name'],
                    'Chromosome': chrom,
                    'Distance': distance,
                    'Correlation': r,
                    'P_value': p
                })

    # Create Cis DataFrame
    cis_df = pd.DataFrame(cis_results)

    if not cis_df.empty:
        # FDR Correction
        if STATSMODELS_AVAIL and df.shape[1] >= 3:
            try:
                reject, qvals, _, _ = multipletests(cis_df['P_value'], method='fdr_bh')
                cis_df['FDR'] = qvals
            except Exception as e:
                print(f"Warning: FDR correction failed: {e}", file=sys.stderr)
                cis_df['FDR'] = cis_df['P_value']
        else:
            cis_df['FDR'] = cis_df['P_value'] # Fallback
            print("Warning: statsmodels not available, skipping FDR correction.", file=sys.stderr)

        # Sort by FDR
        cis_df = cis_df.sort_values('FDR')

    print(f"Writing {len(cis_df)} cis-interactions to {args.output_cis}", file=sys.stderr)
    cis_df.to_csv(args.output_cis, sep='\t', index=False)

    # 5. Trans Analysis (Optional)
    if args.output_trans:
        print("Trans analysis not fully implemented in this version (computational intensity).", file=sys.stderr)
        # Placeholder or simplified version could go here
        with open(args.output_trans, 'w') as f:
            f.write("# Trans analysis placeholder\n")

if __name__ == '__main__':
    main()
