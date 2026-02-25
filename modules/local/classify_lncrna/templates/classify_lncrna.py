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

import os
import sys
from collections import defaultdict
import platform


def parse_gtf(gtf_file):
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

            if feature == 'gene':
                genes[gene_name]['strand'] = strand
                genes[gene_name]['chrom'] = chrom
                genes[gene_name]['start'] = start
                genes[gene_name]['end'] = end

            elif feature == 'exon':
                genes[gene_name]['exons'].append((start, end))
                genes[gene_name]['strand'] = strand
                genes[gene_name]['chrom'] = chrom

            elif feature == 'transcript':
                transcripts[transcript_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_name': gene_name,
                    'attributes': attr_dict
                }

    for gene_name, info in genes.items():
        if ('start' not in info or 'end' not in info) and info['exons']:
            info['start'] = min(e[0] for e in info['exons'])
            info['end'] = max(e[1] for e in info['exons'])

    return genes, transcripts


def find_closest_gene(lncrna_chrom, lncrna_start, lncrna_end, protein_genes):
    min_distance = float('inf')
    closest_gene = None

    for gene_name, gene_info in protein_genes.items():
        if gene_info.get('chrom') != lncrna_chrom:
            continue

        gene_start = gene_info.get('start')
        gene_end = gene_info.get('end')

        if gene_start is None or gene_end is None:
            continue

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


def check_exon_overlap(lncrna_start, lncrna_end, exons):
    for exon_start, exon_end in exons:
        if not (lncrna_end < exon_start or lncrna_start > exon_end):
            return True
    return False


def classify_lncrna(lncrna_info, closest_gene_name, distance, protein_genes):
    lncrna_strand = lncrna_info['strand']
    lncrna_start = lncrna_info['start']
    lncrna_end = lncrna_info['end']

    if closest_gene_name is None:
        return "Intergenic"

    gene_info = protein_genes[closest_gene_name]
    gene_strand = gene_info['strand']

    if distance > 1000:
        return "Intergenic"
    elif 0 < distance <= 1000:
        if gene_strand != lncrna_strand:
            return "Bidirectional"
        else:
            return "Intergenic"
    elif distance == 0:
        if gene_strand != lncrna_strand:
            return "Antisense"
        else:
            exons = gene_info.get('exons', [])
            if check_exon_overlap(lncrna_start, lncrna_end, exons):
                return "Exonic Sense"
            else:
                return "Intronic Sense"

    return "Unknown"


def main():
    lncrna_gtf = "$lncrna_gtf"
    protein_gtf = "$protein_coding_gtf"
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "${meta.id}"
    output_path = f"{prefix}.classification.txt"
    stats_path = f"{prefix}.stats.txt"

    print("Parsing protein-coding genes...", file=sys.stderr)
    protein_genes, _ = parse_gtf(protein_gtf)

    print("Parsing lncRNAs...", file=sys.stderr)
    _, lncrna_transcripts = parse_gtf(lncrna_gtf)

    print("Classifying lncRNAs...", file=sys.stderr)
    classifications = {}
    category_counts = defaultdict(int)

    for transcript_id, lncrna_info in lncrna_transcripts.items():
        closest_gene, distance = find_closest_gene(
            lncrna_info['chrom'],
            lncrna_info['start'],
            lncrna_info['end'],
            protein_genes
        )

        category = classify_lncrna(lncrna_info, closest_gene, distance, protein_genes)

        lncrna_name = lncrna_info.get('gene_name', transcript_id)
        classifications[transcript_id] = {
            'gene_name': lncrna_name,
            'category': category,
            'closest_gene': closest_gene or 'None',
            'distance': distance
        }
        category_counts[category] += 1

    print(f"Writing classification to {output_path}...", file=sys.stderr)
    with open(output_path, 'w') as out:
        out.write("Transcript_ID\tGene_ID\tCategory\tClosest_Gene\tDistance\tChromosome\tStart\tEnd\tStrand\\n")
        for transcript_id, info in sorted(classifications.items()):
            ln = lncrna_transcripts[transcript_id]
            out.write(
                f"{transcript_id}\t{info['gene_name']}\t{info['category']}\t"
                f"{info['closest_gene']}\t{info['distance']}\t"
                f"{ln['chrom']}\t{ln['start']}\t{ln['end']}\t{ln['strand']}\\n"
            )

    print(f"Writing statistics to {stats_path}...", file=sys.stderr)
    total = len(classifications)
    with open(stats_path, 'w') as stats:
        stats.write("lncRNA Classification Summary\\n")
        stats.write("=" * 50 + "\\n\\n")
        stats.write(f"Total lncRNAs: {total}\\n\\n")
        stats.write("Category Breakdown:\\n")
        stats.write("-" * 50 + "\\n")
        for category in sorted(category_counts.keys()):
            count = category_counts[category]
            percentage = (count / total * 100) if total > 0 else 0
            stats.write(f"{category:20s}: {count:6d} ({percentage:5.1f}%)\\n")

    out_dir = os.path.dirname(output_path) or "."
    readme_path = os.path.join(out_dir, "README.txt")
    print(f"Writing README to {readme_path}...", file=sys.stderr)
    with open(readme_path, 'w') as readme:
        readme.write(
            "This file describes the columns of lncrnas_classification.classification.txt\\n\\n"
            "Columns:\\n"
            "- Transcript_ID: lncRNA transcript identifier.\\n"
            "- Gene_ID: lncRNA gene identifier or name.\\n"
            "- Category: genomic classification of the lncRNA relative to the closest protein-coding gene.\\n"
            "- Closest_Gene: name of the closest protein-coding gene.\\n"
            "- Distance: minimal genomic distance in base pairs between the lncRNA and the closest gene.\\n\\n"
            "- Chromosome: chromosome where the lncRNA is located.\\n"
            "- Start: 1-based genomic start position of the lncRNA.\\n"
            "- End: 1-based genomic end position of the lncRNA.\\n"
            "- Strand: genomic strand of the lncRNA ('+' or '-').\\n\\n"
            "Distance rules used for classification:\\n"
            "- > 1000 bp: Intergenic\\n"
            "- 0 < distance ≤ 1000 bp and opposite strand: Bidirectional\\n"
            "- distance = 0 and opposite strand: Antisense\\n"
            "- distance = 0, same strand and overlapping exons: Exonic Sense\\n"
            "- distance = 0, same strand and no exon overlap: Intronic Sense\\n"
        )

    print("Done!", file=sys.stderr)


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
