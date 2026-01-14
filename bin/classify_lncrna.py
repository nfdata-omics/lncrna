#!/usr/bin/env python3

"""
Classify lncRNAs based on their genomic location relative to protein-coding genes.

Categories:
- Intergenic: > 1kb from any protein-coding gene
- Bidirectional: 0-1kb from gene on opposite strand
- Antisense: Overlaps gene on opposite strand
- Exonic Sense: Overlaps exon(s) of gene on same strand
- Intronic Sense: Within introns of gene on same strand
"""

import argparse
import sys
from collections import defaultdict


def parse_gtf(gtf_file):
    """Parse GTF file and extract gene/exon information."""
    genes = defaultdict(lambda: {'strand': None, 'exons': []})
    transcripts = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            start, end = int(start), int(end)

            # Parse attributes
            attr_dict = {}
            for attr in attributes.strip().rstrip(';').split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')

            gene_name = attr_dict.get('gene_name') or attr_dict.get('gene_id')
            transcript_id = attr_dict.get('transcript_id')

            if not gene_name:
                continue

            # Store gene info
            if feature == 'gene':
                genes[gene_name]['strand'] = strand
                genes[gene_name]['chrom'] = chrom
                genes[gene_name]['start'] = start
                genes[gene_name]['end'] = end

            # Store exon info
            elif feature == 'exon':
                genes[gene_name]['exons'].append((start, end))
                genes[gene_name]['strand'] = strand
                genes[gene_name]['chrom'] = chrom

            # Store transcript info for lncRNAs
            elif feature == 'transcript':
                transcripts[transcript_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_name': gene_name,
                    'attributes': attr_dict
                }

    return genes, transcripts


def find_closest_gene(lncrna_chrom, lncrna_start, lncrna_end, protein_genes):
    """Find the closest protein-coding gene to a lncRNA."""
    min_distance = float('inf')
    closest_gene = None

    for gene_name, gene_info in protein_genes.items():
        if gene_info.get('chrom') != lncrna_chrom:
            continue

        gene_start = gene_info.get('start')
        gene_end = gene_info.get('end')

        if gene_start is None or gene_end is None:
            continue

        # Calculate distance
        if lncrna_end < gene_start:
            distance = gene_start - lncrna_end
        elif lncrna_start > gene_end:
            distance = lncrna_start - gene_end
        else:
            distance = 0  # Overlapping

        if distance < min_distance:
            min_distance = distance
            closest_gene = gene_name

    return closest_gene, min_distance


def check_exon_overlap(lncrna_start, lncrna_end, exons):
    """Check if lncRNA overlaps with any exons."""
    for exon_start, exon_end in exons:
        if not (lncrna_end < exon_start or lncrna_start > exon_end):
            return True
    return False


def classify_lncrna(lncrna_info, closest_gene_name, distance, protein_genes):
    """Classify lncRNA based on location relative to closest gene."""
    lncrna_strand = lncrna_info['strand']
    lncrna_start = lncrna_info['start']
    lncrna_end = lncrna_info['end']

    if closest_gene_name is None:
        return "Intergenic"

    gene_info = protein_genes[closest_gene_name]
    gene_strand = gene_info['strand']

    # Rule 1: Intergenic (> 1kb away)
    if distance > 1000:
        return "Intergenic"

    # Rule 2: Bidirectional (0-1kb, opposite strand)
    elif 0 < distance <= 1000:
        if gene_strand != lncrna_strand:
            return "Bidirectional"
        else:
            return "Intergenic"

    # Rule 3-5: Overlapping (distance == 0)
    elif distance == 0:
        # Rule 3: Antisense (opposite strand)
        if gene_strand != lncrna_strand:
            return "Antisense"

        # Same strand - check exon overlap
        else:
            exons = gene_info.get('exons', [])
            if check_exon_overlap(lncrna_start, lncrna_end, exons):
                return "Exonic Sense"
            else:
                return "Intronic Sense"

    return "Unknown"


def main():
    parser = argparse.ArgumentParser(description='Classify lncRNAs based on genomic location')
    parser.add_argument('--lncrna_gtf', required=True, help='lncRNA GTF file')
    parser.add_argument('--protein_gtf', required=True, help='Protein-coding genes GTF file')
    parser.add_argument('--output', required=True, help='Output classification file')
    parser.add_argument('--stats', required=True, help='Output statistics file')
    args = parser.parse_args()

    # Parse protein-coding genes
    print("Parsing protein-coding genes...", file=sys.stderr)
    protein_genes, _ = parse_gtf(args.protein_gtf)

    # Parse lncRNAs
    print("Parsing lncRNAs...", file=sys.stderr)
    _, lncrna_transcripts = parse_gtf(args.lncrna_gtf)

    # Classify each lncRNA
    print("Classifying lncRNAs...", file=sys.stderr)
    classifications = {}
    category_counts = defaultdict(int)

    for transcript_id, lncrna_info in lncrna_transcripts.items():
        # Find closest protein-coding gene
        closest_gene, distance = find_closest_gene(
            lncrna_info['chrom'],
            lncrna_info['start'],
            lncrna_info['end'],
            protein_genes
        )

        # Classify
        category = classify_lncrna(lncrna_info, closest_gene, distance, protein_genes)

        lncrna_name = lncrna_info.get('gene_name', transcript_id)
        classifications[lncrna_name] = {
            'category': category,
            'closest_gene': closest_gene or 'None',
            'distance': distance
        }
        category_counts[category] += 1

    # Write classification output
    print(f"Writing classification to {args.output}...", file=sys.stderr)
    with open(args.output, 'w') as out:
        out.write("lncRNA_ID\tCategory\tClosest_Gene\tDistance\n")
        for lncrna_name, info in sorted(classifications.items()):
            out.write(f"{lncrna_name}\t{info['category']}\t{info['closest_gene']}\t{info['distance']}\n")

    # Write statistics
    print(f"Writing statistics to {args.stats}...", file=sys.stderr)
    total = len(classifications)
    with open(args.stats, 'w') as stats:
        stats.write("lncRNA Classification Summary\n")
        stats.write("=" * 50 + "\n\n")
        stats.write(f"Total lncRNAs: {total}\n\n")
        stats.write("Category Breakdown:\n")
        stats.write("-" * 50 + "\n")
        for category in sorted(category_counts.keys()):
            count = category_counts[category]
            percentage = (count / total * 100) if total > 0 else 0
            stats.write(f"{category:20s}: {count:6d} ({percentage:5.1f}%)\n")

    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
