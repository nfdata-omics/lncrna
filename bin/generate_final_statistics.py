#!/usr/bin/env python3
"""
Generate final statistics and plots comparing lncRNAs and protein-coding genes.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import sys

def parse_gtf(gtf_file):
"""Parse GTF and count genes/transcripts."""
genes = set()
transcripts = set()

with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        feature = fields
        attributes = fields

        # Extract gene_id and transcript_id
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('gene_id'):
                gene_id = attr.split('"')
                genes.add(gene_id)
            elif attr.startswith('transcript_id'):
                transcript_id = attr.split('"')
                transcripts.add(transcript_id)

return len(genes), len(transcripts)
def parse_cpat_results(cpat_file):
"""Parse CPAT results."""
df = pd.read_csv(cpat_file, sep='\t')
return df

def parse_classification(class_file):
"""Parse classification file."""
df = pd.read_csv(class_file, sep='\t')
return df

def create_cpat_comparison_plot(cpat_lncrna_df, cpat_coding_df, output_file):
"""Create boxplot comparing CPAT scores."""
# Prepare data
cpat_lncrna_df['Type'] = 'lncRNA'
cpat_coding_df['Type'] = 'Protein-coding'

combined = pd.concat([
    cpat_lncrna_df[['coding_prob', 'Type']],
    cpat_coding_df[['coding_prob', 'Type']]
])

# Create plot
plt.figure(figsize=(8, 6))
sns.boxplot(data=combined, x='Type', y='coding_prob', palette=['#3498db', '#e74c3c'])
plt.ylabel('CPAT Coding Probability', fontsize=12)
plt.xlabel('')
plt.title('Coding Potential: lncRNAs vs Protein-coding Genes', fontsize=14, fontweight='bold')
plt.ylim(-0.05, 1.05)
plt.axhline(y=0.364, color='red', linestyle='--', label='Threshold (0.364)')
plt.legend()
plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()
def create_classification_barplot(class_df, output_file):
"""Create barplot of lncRNA categories."""
# Count by category
category_counts = class_df['Category'].value_counts()

# Create plot
plt.figure(figsize=(10, 6))
colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6']
ax = category_counts.plot(kind='bar', color=colors[:len(category_counts)])
plt.ylabel('Number of lncRNAs', fontsize=12)
plt.xlabel('Category', fontsize=12)
plt.title('LncRNA Classification by Genomic Location', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right')

# Add value labels on bars
for i, v in enumerate(category_counts):
    ax.text(i, v + max(category_counts)*0.02, str(v), ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()
def generate_html_report(stats_dict, prefix):
"""Generate HTML report."""
html = f"""
<!DOCTYPE html>
<html>
<head>
<title>LncRNA Analysis Final Report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 40px; }}
h1 {{ color: #2c3e50; }}
h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
th {{ background-color: #3498db; color: white; }}
tr:nth-child(even) {{ background-color: #f2f2f2; }}
.metric {{ font-size: 18px; font-weight: bold; color: #2c3e50; }}
.value {{ font-size: 18px; color: #3498db; }}
img {{ max-width: 800px; margin: 20px 0; border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>LncRNA Identification Pipeline - Final Report</h1>

    <h2>Summary Statistics</h2>
    <table>
        <tr>
            <th>Metric</th>
            <th>Value</th>
        </tr>
        <tr>
            <td class="metric">Total Genes (Final Annotation)</td>
            <td class="value">{stats_dict['total_genes']}</td>
        </tr>
        <tr>
            <td class="metric">Total Transcripts (Final Annotation)</td>
            <td class="value">{stats_dict['total_transcripts']}</td>
        </tr>
        <tr>
            <td class="metric">Novel LncRNA Genes</td>
            <td class="value">{stats_dict['lncrna_genes']}</td>
        </tr>
        <tr>
            <td class="metric">Novel LncRNA Transcripts</td>
            <td class="value">{stats_dict['lncrna_transcripts']}</td>
        </tr>
    </table>

    <h2>CPAT Score Comparison</h2>
    <table>
        <tr>
            <th>Category</th>
            <th>Mean Score</th>
            <th>Median Score</th>
            <th>Std Dev</th>
        </tr>
        <tr>
            <td class="metric">LncRNAs</td>
            <td>{stats_dict['cpat_lncrna_mean']:.3f}</td>
            <td>{stats_dict['cpat_lncrna_median']:.3f}</td>
            <td>{stats_dict['cpat_lncrna_std']:.3f}</td>
        </tr>
        <tr>
            <td class="metric">Protein-coding</td>
            <td>{stats_dict['cpat_coding_mean']:.3f}</td>
            <td>{stats_dict['cpat_coding_median']:.3f}</td>
            <td>{stats_dict['cpat_coding_std']:.3f}</td>
        </tr>
    </table>
    <img src="{prefix}.cpat_comparison.png" alt="CPAT Comparison">

    <h2>LncRNA Classification</h2>
    <table>
        <tr>
            <th>Category</th>
            <th>Count</th>
            <th>Percentage</th>
        </tr>
"""

for category, count in stats_dict['classification_counts'].items():
    percentage = (count / stats_dict['lncrna_genes'] * 100) if stats_dict['lncrna_genes'] > 0 else 0
    html += f"""
        <tr>
            <td>{category}</td>
            <td>{count}</td>
            <td>{percentage:.1f}%</td>
        </tr>
"""

html += f"""
    </table>
    <img src="{prefix}.classification_barplot.png" alt="Classification">

    <h2>Pipeline Completion</h2>
    <p>Analysis completed successfully. All files are available in the output directory.</p>
</body>
</html>
"""

with open(f"{prefix}.final_report.html", 'w') as f:
    f.write(html)
def main():
parser = argparse.ArgumentParser(description='Generate final statistics and report')
parser.add_argument('--final_gtf', required=True, help='Final combined GTF')
parser.add_argument('--lncrna_gtf', required=True, help='LncRNA GTF')
parser.add_argument('--cpat_lncrna', required=True, help='CPAT results for lncRNAs')
parser.add_argument('--cpat_coding', required=True, help='CPAT results for protein-coding')
parser.add_argument('--classification', required=True, help='Classification file')
parser.add_argument('--prefix', required=True, help='Output prefix')
args = parser.parse_args()

print("Parsing GTF files...", file=sys.stderr)
total_genes, total_transcripts = parse_gtf(args.final_gtf)
lncrna_genes, lncrna_transcripts = parse_gtf(args.lncrna_gtf)

print("Parsing CPAT results...", file=sys.stderr)
cpat_lncrna_df = parse_cpat_results(args.cpat_lncrna)
cpat_coding_df = parse_cpat_results(args.cpat_coding)

print("Parsing classification...", file=sys.stderr)
class_df = parse_classification(args.classification)
classification_counts = class_df['Category'].value_counts().to_dict()

# Calculate statistics
stats_dict = {
    'total_genes': total_genes,
    'total_transcripts': total_transcripts,
    'lncrna_genes': lncrna_genes,
    'lncrna_transcripts': lncrna_transcripts,
    'cpat_lncrna_mean': cpat_lncrna_df['coding_prob'].mean(),
    'cpat_lncrna_median': cpat_lncrna_df['coding_prob'].median(),
    'cpat_lncrna_std': cpat_lncrna_df['coding_prob'].std(),
    'cpat_coding_mean': cpat_coding_df['coding_prob'].mean(),
    'cpat_coding_median': cpat_coding_df['coding_prob'].median(),
    'cpat_coding_std': cpat_coding_df['coding_prob'].std(),
    'classification_counts': classification_counts
}

# Generate plots
print("Generating plots...", file=sys.stderr)
create_cpat_comparison_plot(cpat_lncrna_df, cpat_coding_df, f"{args.prefix}.cpat_comparison.png")
create_classification_barplot(class_df, f"{args.prefix}.classification_barplot.png")

# Generate HTML report
print("Generating HTML report...", file=sys.stderr)
generate_html_report(stats_dict, args.prefix)

# Generate TSV summary
print("Writing summary TSV...", file=sys.stderr)
with open(f"{args.prefix}.statistics_summary.tsv", 'w') as f:
    f.write("Metric\tValue\n")
    for key, value in stats_dict.items():
        if key != 'classification_counts':
            f.write(f"{key}\t{value}\n")

print("Done!", file=sys.stderr)
if name == 'main':
main()
