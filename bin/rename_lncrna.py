#!/usr/bin/env python3
"""
Rename lncRNAs based on their closest protein-coding gene.
Example: MSTRG.1.1 → LINC-BRCA1-001
"""

import argparse
import sys
from collections import defaultdict


def parse_gtf_genes(gtf_file):
    """Parse GTF and extract gene locations."""
    genes = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature != 'gene':
                continue

            # Parse attributes
            attr_dict = {}
            for attr in attributes.strip().rstrip(';').split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')

            gene_name = attr_dict.get('gene_name') or attr_dict.get('gene_id')
            if gene_name:
                genes[gene_name] = {
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand
                }

    return genes


def find_closest_gene(lncrna_chrom, lncrna_start, lncrna_end, protein_genes):
    """Find closest protein-coding gene."""
    min_distance = float('inf')
    closest_gene = None

    for gene_name, gene_info in protein_genes.items():
        if gene_info['chrom'] != lncrna_chrom:
            continue

        gene_start = gene_info['start']
        gene_end = gene_info['end']

        # Calculate distance
        if lncrna_end < gene_start:
            distance = gene_start - lncrna_end
        elif lncrna_start > gene_end:
            distance = lncrna_start - gene_end
        else:
            distance = 0

        if distance < min_distance:
            min_distance = distance
            closest_gene = gene_name

    return closest_gene, min_distance


def rename_lncrnas(lncrna_gtf, protein_gtf, output_gtf, mapping_file, prefix):
    """Rename lncRNAs based on closest gene and write mapping."""

    # Parse protein-coding genes
    print("Parsing protein-coding genes...", file=sys.stderr)
    protein_genes = parse_gtf_genes(protein_gtf)

    # Track lncRNA genes and assign new names
    lncrna_genes = {}
    gene_counter = defaultdict(int)  # Counter per protein gene
    mapping = {}

    print("Processing lncRNAs...", file=sys.stderr)

    # First pass: identify unique lncRNA genes and assign names
    with open(lncrna_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end = fields[0], fields[1], fields[2], fields[3], fields[4]
            strand, attributes = fields[6], fields[8]

            # Parse attributes
            attr_dict = {}
            for attr in attributes.strip().rstrip(';').split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')

            gene_id = attr_dict.get('gene_id')

            if not gene_id or gene_id in lncrna_genes:
                continue

            # Find closest protein-coding gene
            closest_gene, distance = find_closest_gene(
                chrom, int(start), int(end), protein_genes
            )

            if closest_gene:
                # Generate new name
                gene_counter[closest_gene] += 1
                new_name = f"{prefix}-{closest_gene}-{gene_counter[closest_gene]:03d}"
            else:
                # No nearby gene, use chromosome
                gene_counter[f"CHR{chrom}"] += 1
                new_name = f"{prefix}-CHR{chrom}-{gene_counter[f'CHR{chrom}']:03d}"

            lncrna_genes[gene_id] = new_name
            mapping[gene_id] = {
                'new_name': new_name,
                'closest_gene': closest_gene or 'None',
                'distance': distance
            }

    # Second pass: write renamed GTF
    print(f"Writing renamed GTF to {output_gtf}...", file=sys.stderr)
    with open(lncrna_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                outfile.write(line)
                continue

            attributes = fields[8]

            # Replace gene_id and gene_name
            for old_id, new_name in lncrna_genes.items():
                if f'gene_id "{old_id}"' in attributes:
                    attributes = attributes.replace(f'gene_id "{old_id}"', f'gene_id "{new_name}"')
                    attributes = attributes.replace(f'gene_name "{old_id}"', f'gene_name "{new_name}"')

                    # Also update transcript_id if present
                    if 'transcript_id' in attributes:
                        import re
                        transcript_pattern = re.compile(r'transcript_id "([^"]+)"')
                        match = transcript_pattern.search(attributes)
                        if match:
                            old_transcript = match.group(1)
                            # Keep transcript suffix if exists
                            if '.' in old_transcript:
                                suffix = old_transcript.split('.')[-1]
                                new_transcript = f"{new_name}.{suffix}"
                            else:
                                new_transcript = new_name
                            attributes = attributes.replace(
                                f'transcript_id "{old_transcript}"',
                                f'transcript_id "{new_transcript}"'
                            )

            fields[8] = attributes
            outfile.write('\t'.join(fields) + '\n')

    # Write mapping file
    print(f"Writing mapping to {mapping_file}...", file=sys.stderr)
    with open(mapping_file, 'w') as mapfile:
        mapfile.write("Old_ID\tNew_Name\tClosest_Gene\tDistance\n")
        for old_id, info in sorted(mapping.items()):
            mapfile.write(
                f"{old_id}\t{info['new_name']}\t{info['closest_gene']}\t{info['distance']}\n"
            )

    print(f"Renamed {len(lncrna_genes)} lncRNA genes.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Rename lncRNAs based on closest gene')
    parser.add_argument('--lncrna_gtf', required=True, help='Input lncRNA GTF')
    parser.add_argument('--protein_gtf', required=True, help='Protein-coding genes GTF')
    parser.add_argument('--output', required=True, help='Output renamed GTF')
    parser.add_argument('--mapping', required=True, help='Output mapping file')
    parser.add_argument('--prefix', default='LINC', help='Prefix for new names (default: LINC)')
    args = parser.parse_args()

    rename_lncrnas(
        args.lncrna_gtf,
        args.protein_gtf,
        args.output,
        args.mapping,
        args.prefix
    )


if __name__ == '__main__':
    main()
