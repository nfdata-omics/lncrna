#!/usr/bin/env python3
"""
Generate comprehensive HTML report for lncRNA pipeline results.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path
import base64
import json


def encode_image(image_path):
    """Encode image to base64 for embedding in HTML."""
    if not Path(image_path).exists():
        return None
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def parse_de_results(de_file):
    """Parse differential expression results."""
    if not de_file or de_file == 'NO_FILE' or not Path(de_file).exists():
        return None, None
    df = pd.read_csv(de_file, sep='\t', index_col=0)
    stats = {
        'total': len(df),
        'significant': sum((df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1)),
        'upregulated': sum((df['padj'] < 0.05) & (df['log2FoldChange'] > 1)),
        'downregulated': sum((df['padj'] < 0.05) & (df['log2FoldChange'] < -1)),
        'lncrna_de': sum((df['padj'] < 0.05) & (df['gene_type'] == 'lncRNA')) if 'gene_type' in df.columns else 0
    }
    return df, stats


def parse_gtf_stats(gtf_file):
    """Parse GTF file to get basic statistics."""
    if not gtf_file or gtf_file == 'NO_FILE' or not Path(gtf_file).exists():
        return {}
    # Simple GTF parsing for basic stats
    gene_types = {}
    total_genes = 0
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                total_genes += 1
                # Extract gene_type from attributes
                attrs = parts[8]
                if 'gene_type' in attrs:
                    gene_type = attrs.split('gene_type "')[1].split('"')[0]
                    gene_types[gene_type] = gene_types.get(gene_type, 0) + 1
    return {'total_genes': total_genes, 'gene_types': gene_types}

def generate_html_report(args):
    """Generate comprehensive HTML report."""
    # Parse input files
    print("Parsing differential expression results...", file=sys.stderr)
    de_df, de_stats = parse_de_results(args.de_results)
    print("Parsing classification...", file=sys.stderr)
    class_df = pd.read_csv(args.classification, sep='\t')
    class_stats_df = pd.read_csv(args.classification_stats, sep='\t')
    print("Parsing GTF statistics...", file=sys.stderr)
    gtf_stats = parse_gtf_stats(args.final_gtf)
    # Parse count summary if available
    count_summary = None
    if args.count_summary and args.count_summary != 'NO_FILE' and Path(args.count_summary).exists():
        count_summary = pd.read_csv(args.count_summary, sep='\t')
    # Parse rename mapping if available
    rename_mapping = None
    if args.rename_mapping and args.rename_mapping != 'NO_FILE' and Path(args.rename_mapping).exists():
        rename_mapping = pd.read_csv(args.rename_mapping, sep='\t')
    # Embed plots
    volcano_img = None
    pca_img = None
    heatmap_img = None
    ma_img = None
    if args.de_plots and args.de_plots != 'NO_FILE':
        de_plots_path = Path(args.de_plots)
        if de_plots_path.is_dir():
            # Handle directory of plots
            volcano_img = encode_image(de_plots_path / "volcano_plot.png")
            pca_img = encode_image(de_plots_path / "pca_plot.png")
            heatmap_img = encode_image(de_plots_path / "heatmap.png")
            ma_img = encode_image(de_plots_path / "ma_plot.png")
        else:
            # Handle single plot files (you might need to adjust this)
            volcano_img = encode_image(args.de_plots) if 'volcano' in str(args.de_plots) else None
    cpat_img = encode_image(args.cpat_comparison) if args.cpat_comparison != 'NO_FILE' else None
    # Generate DE section HTML
    de_section = ""
    if de_df is not None and de_stats is not None:
        de_section = f"""
        <h2> Differential Expression Analysis</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-label">Total Genes</div>
                <div class="stat-value">{de_stats['total']:,}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">DE Genes</div>
                <div class="stat-value">{de_stats['significant']:,}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Upregulated</div>
                <div class="stat-value">{de_stats['upregulated']:,}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Downregulated</div>
                <div class="stat-value">{de_stats['downregulated']:,}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">DE LncRNAs</div>
                <div class="stat-value">{de_stats['lncrna_de']:,}</div>
            </div>
        </div>
        """
        # Add DE plots if available
        if volcano_img or pca_img:
            de_section += '<div class="grid-2">'
            if volcano_img:
                de_section += f'''
                <div class="plot-container">
                    <h3>Volcano Plot</h3>
                    <img src="data:image/png;base64,{volcano_img}" alt="Volcano Plot">
                </div>
                '''
            if pca_img:
                de_section += f'''
                <div class="plot-container">
                    <h3>PCA Plot</h3>
                    <img src="data:image/png;base64,{pca_img}" alt="PCA Plot">
                </div>
                '''
            de_section += '</div>'
        if heatmap_img:
            de_section += f'''
            <div class="plot-container">
                <h3>Heatmap - Top DE Genes</h3>
                <img src="data:image/png;base64,{heatmap_img}" alt="Heatmap">
            </div>
            '''
        # Top DE lncRNAs table
        if 'gene_type' in de_df.columns:
            lncrna_de = de_df[de_df['gene_type'] == 'lncRNA'].head(20)
            if not lncrna_de.empty:
                de_section += f'''
                <h3>Top Differentially Expressed LncRNAs</h3>
                {lncrna_de.to_html(classes='table', index=True)}
                '''
    else:
        de_section = '''
        <h2> Differential Expression Analysis</h2>
        <p><em>Differential expression analysis was not performed or no design file was provided.</em></p>
        '''
    # Generate pipeline summary stats
    pipeline_stats = f"""
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-label">Total Genes in Final GTF</div>
            <div class="stat-value">{gtf_stats.get('total_genes', 'N/A'):,}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">LncRNA Classifications</div>
            <div class="stat-value">{len(class_df):,}</div>
        </div>
    """
    if rename_mapping is not None:
        pipeline_stats += f'''
        <div class="stat-card">
            <div class="stat-label">Renamed LncRNAs</div>
            <div class="stat-value">{len(rename_mapping):,}</div>
        </div>
        '''
    pipeline_stats += "</div>"
    # Generate HTML
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>LncRNA Pipeline - Final Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 0;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: white;
        }}
        header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        h1 {{ margin: 0; font-size: 2.5em; }}
        h2 {{ color: #667eea; border-bottom: 3px solid #667eea; padding-bottom: 10px; margin-top: 40px; }}
        h3 {{ color: #764ba2; margin-top: 30px; }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .stat-value {{
            font-size: 3em;
            font-weight: bold;
            margin: 10px 0;
        }}
        .stat-label {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #667eea;
            color: white;
            font-weight: bold;
        }}
        tr:hover {{ background-color: #f5f5f5; }}
        .plot-container {{
            margin: 30px 0;
            text-align: center;
        }}
        .plot-container img {{
            max-width: 100%;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .grid-2 {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 30px;
            margin: 30px 0;
        }}
        footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            margin-top: 50px;
            border-top: 1px solid #ddd;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1> LncRNA Identification Pipeline</h1>
            <p style="font-size: 1.2em; margin-top: 10px;">Comprehensive Analysis Report</p>
        </header>

        <h2> Pipeline Summary</h2>
        {pipeline_stats}
        {de_section}
        <h2> LncRNA Classification</h2>
        <h3>Classification Statistics</h3>
        {class_stats_df.to_html(classes='table', index=False)}
        <h3>Individual Classifications</h3>
        {class_df.head(50).to_html(classes='table', index=False)}
    """
    if cpat_img:
        html += f'''
        <div class="plot-container">
            <h3>CPAT Score Comparison</h3>
            <img src="data:image/png;base64,{cpat_img}" alt="CPAT Comparison">
        </div>
        '''
    if count_summary is not None:
        html += f'''
        <h2> Sample Statistics</h2>
        {count_summary.to_html(classes='table', index=False)}
        '''
    if rename_mapping is not None:
        html += f'''
        <h2> LncRNA Renaming</h2>
        <h3>Rename Mapping (First 50)</h3>
        {rename_mapping.head(50).to_html(classes='table', index=False)}
        '''
    html += '''
        <footer>
            <p>Generated by nf-core/lncrna pipeline</p>
            <p>For questions or issues, please visit the GitHub repository</p>
        </footer>
    </div>
</body>
</html>
'''

    # Write HTML report
    with open(args.output, 'w') as f:
        f.write(html)
    # Write text summary
    with open(args.summary, 'w') as f:
        f.write("=== LncRNA Pipeline Summary ===\n\n")
        f.write(f"Total genes in final GTF: {gtf_stats.get('total_genes', 'N/A')}\n")
        f.write(f"LncRNA classifications: {len(class_df)}\n\n")
        if de_stats:
            f.write("Differential Expression Results:\n")
            f.write(f"  Total genes analyzed: {de_stats['total']}\n")
            f.write(f"  Differentially expressed genes: {de_stats['significant']}\n")
            f.write(f"    - Upregulated: {de_stats['upregulated']}\n")
            f.write(f"    - Downregulated: {de_stats['downregulated']}\n")
            f.write(f"  Differentially expressed lncRNAs: {de_stats['lncrna_de']}\n\n")
        else:
            f.write("Differential expression analysis: Not performed\n\n")
        f.write("LncRNA Classification:\n")
        for _, row in class_stats_df.iterrows():
            f.write(f"  - {row.iloc[0]}: {row.iloc[1]}\n")
    # Write JSON stats
    stats_json = {
        'pipeline': 'lncrna',
        'total_genes': gtf_stats.get('total_genes', 0),
        'lncrna_classifications': len(class_df),
        'differential_expression': de_stats if de_stats else {},
        'classification_stats': class_stats_df.to_dict('records')
    }
    json_output = args.output.replace('.html', '.stats.json')
    with open(json_output, 'w') as f:
        json.dump(stats_json, f, indent=2)
    print(f"Report generated: {args.output}", file=sys.stderr)
    print(f"Summary generated: {args.summary}", file=sys.stderr)
    print(f"Stats JSON generated: {json_output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Generate lncRNA pipeline report')
    parser.add_argument('--de_results', help='Differential expression results')
    parser.add_argument('--de_plots', help='Directory with DE plots')
    parser.add_argument('--classification', required=True, help='LncRNA classification file')
    parser.add_argument('--classification_stats', required=True, help='Classification statistics')
    parser.add_argument('--cpat_comparison', required=True, help='CPAT comparison plot')
    parser.add_argument('--count_summary', help='Count summary file')
    parser.add_argument('--final_gtf', required=True, help='Final GTF annotation')
    parser.add_argument('--rename_mapping', help='Rename mapping file')
    parser.add_argument('--output', required=True, help='Output HTML report')
    parser.add_argument('--summary', required=True, help='Output text summary')
    args = parser.parse_args()
    generate_html_report(args)

if __name__ == '__main__':
    main()
