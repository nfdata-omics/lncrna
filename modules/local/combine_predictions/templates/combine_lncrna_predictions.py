#!/usr/bin/env python3

"""
Combine predictions from CPAT, FEELnc, and PLEK
Use consensus voting to identify high-confidence lncRNAs
"""

import sys
from collections import defaultdict
from Bio import SeqIO
import platform


def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def parse_cpat(cpat_file):
    """Parse CPAT results"""
    predictions = {}

    with open(cpat_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                transcript_id = fields[0]
                coding_prob = float(fields[10])
                # CPAT: coding_prob < 0.364 = Non-coding (for human)
                pred = 'Non-coding' if coding_prob < 0.364 else 'Coding'
                predictions[transcript_id] = {
                    'prediction': pred,
                    'score': coding_prob
                }

    return predictions


def parse_feelnc(feelnc_file):
    """Parse FEELnc results"""
    predictions = {}

    with open(feelnc_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                transcript_id = fields[0]
                pred = fields[1]  # 'lncRNA' or 'mRNA'
                score = float(fields[2]) if fields[2] != 'NA' else 0.5

                # Standardize
                pred_class = 'Non-coding' if pred == 'lncRNA' else 'Coding'
                predictions[transcript_id] = {
                    'prediction': pred_class,
                    'score': score
                }

    return predictions


def parse_plek(plek_file):
    """Parse PLEK results"""
    predictions = {}

    with open(plek_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                transcript_id = fields[0]
                pred = fields[1]
                score = fields[2]

                predictions[transcript_id] = {
                    'prediction': pred,
                    'score': score
                }

    return predictions


def consensus_vote(cpat_pred, feelnc_pred, plek_pred, mode='majority'):
    """
    Combine predictions using consensus voting

    Modes:
    - majority: At least 2/3 tools agree on Non-coding
    - strict: All 3 tools agree on Non-coding
    - lenient: At least 1/3 tools predict Non-coding
    """

    votes = []
    if cpat_pred:
        votes.append(cpat_pred['prediction'])
    if feelnc_pred:
        votes.append(feelnc_pred['prediction'])
    if plek_pred:
        votes.append(plek_pred['prediction'])

    non_coding_count = votes.count('Non-coding')
    total_votes = len(votes)

    if mode == 'strict':
        return 'lncRNA' if non_coding_count == total_votes else 'Coding Potential Evidence'
    elif mode == 'majority':
        return 'lncRNA' if non_coding_count >= (total_votes / 2) else 'Coding Potential Evidence'
    elif mode == 'lenient':
        return 'lncRNA' if non_coding_count >= 1 else 'Coding Potential Evidence'
    else:
        return 'lncRNA' if non_coding_count >= 2 else 'Coding Potential Evidence'


def combine_predictions(cpat_file, feelnc_file, plek_file, gtf_file,
                        fasta_file, tmap_file, mode, output_prefix):
    """Combine all predictions and generate final lncRNA set"""

    print("Parsing predictions...", file=sys.stderr)
    cpat_preds = parse_cpat(cpat_file)
    feelnc_preds = parse_feelnc(feelnc_file)
    plek_preds = parse_plek(plek_file)

    # Get all transcript IDs
    all_transcripts = set(cpat_preds.keys()) | set(feelnc_preds.keys()) | set(plek_preds.keys())

    print(f"Total transcripts: {len(all_transcripts)}", file=sys.stderr)

    # Consensus voting
    final_predictions = {}
    stats = {'lncRNA': 0, 'Coding Potential Evidence': 0}

    for transcript_id in all_transcripts:
        cpat = cpat_preds.get(transcript_id)
        feelnc = feelnc_preds.get(transcript_id)
        plek = plek_preds.get(transcript_id)

        consensus = consensus_vote(cpat, feelnc, plek, mode)

        final_predictions[transcript_id] = {
            'consensus': consensus,
            'cpat': cpat['prediction'] if cpat else 'NA',
            'feelnc': feelnc['prediction'] if feelnc else 'NA',
            'plek': plek['prediction'] if plek else 'NA',
            'cpat_score': cpat['score'] if cpat else 'NA',
            'feelnc_score': feelnc['score'] if feelnc else 'NA',
            'plek_score': plek['score'] if plek else 'NA'
        }

        stats[consensus] += 1

    print(f"Consensus results: {stats['lncRNA']} lncRNAs, {stats['Coding Potential Evidence']} Coding Potential Evidence",
          file=sys.stderr)

    # Filter lncRNAs
    lncrna_ids = {tid for tid, pred in final_predictions.items()
                  if pred['consensus'] == 'lncRNA'}

    # Write summary
    with open(f"{output_prefix}.prediction_summary.tsv", 'w') as out:
        out.write("transcript_id\tconsensus\tcpat\tfeelnc\tplek\t"
                  "cpat_score\tfeelnc_score\tplek_score\\n")
        for tid in sorted(all_transcripts):
            pred = final_predictions[tid]
            out.write(f"{tid}\t{pred['consensus']}\t{pred['cpat']}\t"
                      f"{pred['feelnc']}\t{pred['plek']}\t"
                      f"{pred['cpat_score']}\t{pred['feelnc_score']}\t"
                      f"{pred['plek_score']}\\n")

    # Write lncRNA-only summary
    with open(f"{output_prefix}.prediction_summary_lncrna.tsv", 'w') as out:
        out.write("transcript_id\tconsensus\tcpat\tfeelnc\tplek\t"
                  "cpat_score\tfeelnc_score\tplek_score\\n")
        for tid in sorted(lncrna_ids):
            pred = final_predictions[tid]
            out.write(f"{tid}\t{pred['consensus']}\t{pred['cpat']}\t"
                      f"{pred['feelnc']}\t{pred['plek']}\t"
                      f"{pred['cpat_score']}\t{pred['feelnc_score']}\t"
                      f"{pred['plek_score']}\\n")

    # Filter GTF
    print("Filtering GTF...", file=sys.stderr)
    with open(gtf_file, 'r') as gin, \
         open(f"{output_prefix}.final_lncrna.gtf", 'w') as gout:
        for line in gin:
            if line.startswith('#'):
                gout.write(line)
                continue
            if any(f'transcript_id "{tid}"' in line for tid in lncrna_ids):
                gout.write(line)

    # Filter FASTA
    print("Filtering FASTA...", file=sys.stderr)
    with open(f"{output_prefix}.final_lncrna.fa", 'w') as fout:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in lncrna_ids:
                SeqIO.write(record, fout, "fasta")

    # Write report
    with open(f"{output_prefix}.lncrna_report.txt", 'w') as report:
        report.write("=" * 70 + "\\n")
        report.write("LNCRNA IDENTIFICATION REPORT\\n")
        report.write("=" * 70 + "\\n\\n")
        report.write(f"Consensus mode: {mode}\\n\\n")
        report.write(f"Total transcripts analyzed: {len(all_transcripts)}\\n")
        report.write(f"Predicted lncRNAs: {stats['lncRNA']} ({stats['lncRNA']/len(all_transcripts)*100:.2f}%)\\n")
        report.write(f"Predicted Coding Potential Evidence: {stats['Coding Potential Evidence']} ({stats['Coding Potential Evidence']/len(all_transcripts)*100:.2f}%)\\n")

    print(f"Final lncRNA set: {len(lncrna_ids)} transcripts", file=sys.stderr)


def main():
    cpat_file = "$cpat_results"
    feelnc_file = "$feelnc_results"
    plek_file = "$plek_results"
    gtf_file = "$gtf"
    fasta_file = "$fasta"
    tmap_file = "$tmap"

    mode_value = "$task.ext.args" if "$task.ext.args" != "null" else "majority"
    if mode_value not in ["strict", "majority", "lenient"]:
        mode_value = "majority"

    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "${meta4.id}"

    combine_predictions(cpat_file, feelnc_file, plek_file, gtf_file,
                        fasta_file, tmap_file, mode_value, prefix)


if __name__ == '__main__':
    main()
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
