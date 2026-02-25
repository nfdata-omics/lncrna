#!/usr/bin/env python

# Written by Olga Botvinnik with subsequent reworking by Jonathan Manning. Released under the MIT license.

import logging
import re
import statistics
import os
import platform
from typing import Set

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("fasta_gtf_filter")
logger.setLevel(logging.INFO)


def extract_fasta_seq_names(fasta_name: str) -> Set[str]:
    """Extracts the sequence names from a FASTA file."""
    with open(fasta_name) as fasta:
        return {line[1:].split(None, 1)[0] for line in fasta if line.startswith(">")}


def tab_delimited(file: str) -> float:
    """Check if file is tab-delimited and return median number of tabs."""
    with open(file, "r") as f:
        data = f.read(102400)
        return statistics.median(line.count("\\t") for line in data.split("\\n"))


def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def filter_gtf(fasta: str, gtf_in: str, filtered_gtf_out: str, skip_transcript_id_check: bool) -> None:
    """Filter GTF file based on FASTA sequence names."""
    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 9 tab-separated columns.")

    seq_names_in_genome = extract_fasta_seq_names(fasta)
    logger.info(f"Extracted chromosome sequence names from {fasta}")
    logger.debug("All sequence IDs from FASTA: " + ", ".join(sorted(seq_names_in_genome)))

    seq_names_in_gtf = set()
    try:
        with open(gtf_in) as gtf, open(filtered_gtf_out, "w") as out:
            line_count = 0
            for line in gtf:
                seq_name = line.split("\\t")[0]
                seq_names_in_gtf.add(seq_name)  # Add sequence name to the set

                if seq_name in seq_names_in_genome:
                    if skip_transcript_id_check or re.search(r'transcript_id "([^"]+)"', line):
                        out.write(line)
                        line_count += 1

            if line_count == 0:
                raise ValueError("All GTF lines removed by filters")

    except IOError as e:
        logger.error(f"File operation failed: {e}")
        return

    logger.debug("All sequence IDs from GTF: " + ", ".join(sorted(seq_names_in_gtf)))
    logger.info(f"Extracted {line_count} matching sequences from {gtf_in} into {filtered_gtf_out}")


if __name__ == "__main__":
    fasta = "$fasta"
    gtf_in = "$gtf"
    prefix = os.path.splitext(os.path.basename(fasta))[0]
    skip_transcript_id_check = "$task.ext.args" == "--skip_transcript_id_check"
    filtered_gtf_out = prefix + ".filtered.gtf"
    filter_gtf(fasta, gtf_in, filtered_gtf_out, skip_transcript_id_check)
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
