#!/usr/bin/env python3
"""
Parse FEELnc_codpot output, standardize results and split FASTA.
"""

import sys
import platform
from Bio import SeqIO


def parse_feelnc_output(input_file, output_file, fasta_file, prefix):

    results = {}
    print(f"Parsing FEELnc output: {input_file}", file=sys.stderr)

    with open(input_file, 'r') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if header is None:
                header = line.split('\t')
                continue                         # salta el header

            fields = line.split('\t')
            if len(fields) < 11:
                continue

            transcript_id      = fields[0]
            coding_potential   = float(fields[9])
            label              = int(fields[10])

            # label: 1 = coding (mRNA), 0 = non-coding (lncRNA)
            classification = 'Coding Potential Evidence' if label == 1 else 'lncRNA'

            results[transcript_id] = {
                'coding_potential' : coding_potential,
                'classification'   : classification
            }

    print(f"Parsed {len(results)} transcripts", file=sys.stderr)

    # Write standardized TSV
    with open(output_file, 'w') as out:
        out.write("transcript_id\tclassification\tcoding_potential\\n")
        for transcript_id, data in results.items():
            out.write(
                f"{transcript_id}\t{data['classification']}\t"
                f"{data['coding_potential']}\\n"
            )

    # Split into lncRNA / mRNA FASTAs
    lncrna_ids = {tid for tid, d in results.items() if d['classification'] == 'lncRNA'}
    mrna_ids   = {tid for tid, d in results.items() if d['classification'] == 'Coding Potential Evidence'}

    print(f"lncRNAs: {len(lncrna_ids)}, Coding Potential Evidence: {len(mrna_ids)}", file=sys.stderr)

    if lncrna_ids:
        with open(f"{prefix}.lncRNA.fa", 'w') as lnc_out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in lncrna_ids:
                    SeqIO.write(record, lnc_out, "fasta")

    if mrna_ids:
        with open(f"{prefix}.mRNA.fa", 'w') as mrna_out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in mrna_ids:
                    SeqIO.write(record, mrna_out, "fasta")


def main():
    # Variables from Nextflow — NO genome_fasta, NO FEELnc execution
    input_file     = "$feelnc_rf_txt"       # *_RF.txt from FEELNC_CODPOT_RUN
    candidate_fasta = "$candidate_fasta"    # .noORF.fa from FEELNC_CODPOT_RUN
    prefix         = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "${meta.id}"
    output_file    = f"{prefix}.feelnc.tsv"

    try:
        parse_feelnc_output(input_file, output_file, candidate_fasta, prefix)
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        sys.exit(1)

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "feelnc": "0.2.1"
        }
    }
    with open("versions.yml", "w") as f:
        for k, v in versions.items():
            f.write(f"{k}:\\n")
            for k2, v2 in v.items():
                f.write(f"  {k2}: {v2}\\n")


if __name__ == '__main__':
    main()
