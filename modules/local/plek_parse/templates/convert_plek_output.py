#!/usr/bin/env python3

"""
Convert PLEK output to standardized TSV format
"""

import sys
import platform


def format_yaml_like(data: dict, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        if isinstance(value, dict):
            yaml_str += " " * indent + f"{key}:\\n"
            yaml_str += format_yaml_like(value, indent + 4)
        else:
            yaml_str += " " * indent + f"{key}: {value}\\n"
    return yaml_str


def convert_plek_output(input_file, output_file):
    results = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                continue

            classification_raw = fields[0]   # "Non-coding" o "Coding"
            score              = fields[1]   # "-0.687237"
            transcript_id      = fields[2].lstrip('>')  # ">MSTRG.21.1" → "MSTRG.21.1"

            pred_class = 'Non-coding' if 'Non-coding' in classification_raw else 'Coding'

            results.append({
                'transcript_id' : transcript_id,
                'classification': pred_class,
                'score'         : score
            })

    with open(output_file, 'w') as out:
        out.write("transcript_id\tclassification\tscore\\n")
        for result in results:
            out.write(f"{result['transcript_id']}\t{result['classification']}\t{result['score']}\\n")

    print(f"Converted {len(results)} PLEK predictions to {output_file}", file=sys.stderr)

prefix = "${prefix}"
convert_plek_output("${plek_raw}", "${prefix}.plek.tsv")

versions_this_module = {}
versions_this_module["${task.process}"] = {"python": platform.python_version()}
with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions_this_module))
