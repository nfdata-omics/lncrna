#!/usr/bin/env python3

"""
Parse FEELnc_codpot output, standardize results and split FASTA.
"""

import sys
import subprocess
import shlex
import platform
from Bio import SeqIO


def parse_feelnc_output(input_file, output_file, fasta_file, prefix):
    """
    Parse FEELnc output and create standardized TSV

    FEELnc output format:
    transcript_id   lncRNA_probability   mRNA_probability   classification
    """

    results = {}

    print(f"Parsing FEELnc output: {input_file}", file=sys.stderr)

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('transcript'):
                continue

            fields = line.strip().split('\t')
            if len(fields) >= 4:
                transcript_id = fields[0]
                lncrna_prob = float(fields[1])
                mrna_prob = float(fields[2])
                classification = fields[3]

                results[transcript_id] = {
                    'lncrna_prob': lncrna_prob,
                    'mrna_prob': mrna_prob,
                    'classification': classification
                }

    print(f"Parsed {len(results)} transcripts", file=sys.stderr)

    # Write standardized TSV
    with open(output_file, 'w') as out:
        out.write("transcript_id\tclassification\tscore\tlncrna_prob\tmrna_prob\n")

        for transcript_id, data in results.items():
            out.write(f"{transcript_id}\t{data['classification']}\t"
                      f"{data['lncrna_prob']}\t{data['lncrna_prob']}\t"
                      f"{data['mrna_prob']}\n")

    # Separate lncRNA and mRNA sequences
    lncrna_ids = {tid for tid, data in results.items()
                  if data['classification'] == 'lncRNA'}
    mrna_ids = {tid for tid, data in results.items()
                if data['classification'] == 'mRNA'}

    print(f"lncRNAs: {len(lncrna_ids)}, mRNAs: {len(mrna_ids)}", file=sys.stderr)

    # Write lncRNA FASTA
    if lncrna_ids:
        with open(f"{prefix}.lncRNA.fa", 'w') as lnc_out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in lncrna_ids:
                    SeqIO.write(record, lnc_out, "fasta")

    # Write mRNA FASTA
    if mrna_ids:
        with open(f"{prefix}.mRNA.fa", 'w') as mrna_out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in mrna_ids:
                    SeqIO.write(record, mrna_out, "fasta")

    print(f"✓ Created {output_file}", file=sys.stderr)


def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def main():
    candidate_fasta = "$candidate_fasta"
    genome_fasta = "$genome_fasta"
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "${meta.id}"

    mrna_ref = "$task.ext.mrna_ref"
    if not mrna_ref or mrna_ref == "null":
        mrna_ref = ""

    extra_args = "$task.ext.args"
    extra_args_list = shlex.split(extra_args) if extra_args and extra_args != "null" else []

    cmd = [
        "FEELnc_codpot.pl",
        "-i",
        candidate_fasta,
        "-a",
        genome_fasta,
        "--mode=shuffle",
        "--numtx=500",
        "-o",
        f"{prefix}.feelnc_codpot.txt",
    ]

    if mrna_ref:
        cmd.extend(["-m", mrna_ref])

    cmd.extend(extra_args_list)

    print("Running FEELnc_codpot:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

    input_file = f"{prefix}.feelnc_codpot.txt"
    output_file = f"{prefix}.feelnc.tsv"

    try:
        parse_feelnc_output(input_file, output_file, candidate_fasta, prefix)
    except Exception as e:
        print(f"✗ Error parsing FEELnc output: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        feelnc_version = subprocess.check_output(
            "FEELnc_codpot.pl --version 2>&1 | grep -oP 'FEELnc version \\K[0-9.]+' || echo '0.2.1'",
            shell=True,
            text=True,
        ).strip()
    except Exception:
        feelnc_version = "0.2.1"

    versions_this_module = {
        "${task.process}": {
            "python": platform.python_version(),
            "feelnc": feelnc_version,
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))


if __name__ == '__main__':
    main()
