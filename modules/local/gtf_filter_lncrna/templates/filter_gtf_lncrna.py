#!/usr/bin/env python3

"""
Filter GTF file to extract only lncRNA entries.
Designed to work after filter_gtf.py for validation.

Written by Karla Ruiz. Released under the MIT license.
"""

import logging
import argparse
import re
import platform
from typing import Optional, Set

logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("filter_gtf_lncrna")
logger.setLevel(logging.INFO)

# Compile regex once at module level
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_ID_RE       = re.compile(r'gene_id "([^"]+)"')
ATTR_RE          = {}  # cache for dynamic attribute patterns


def extract_attribute(attributes_str: str, attribute_name: str) -> Optional[str]:
    """Extract attribute value from GTF attributes string."""
    if attribute_name not in ATTR_RE:
        ATTR_RE[attribute_name] = re.compile(f'{attribute_name} "([^"]+)"')
    match = ATTR_RE[attribute_name].search(attributes_str)
    return match.group(1) if match else None


def is_lncrna_biotype(attributes_str: str, biotypes: set) -> bool:
    return any(
        f'gene_type "{bt}"'          in attributes_str or
        f'gene_biotype "{bt}"'       in attributes_str or
        f'transcript_type "{bt}"'    in attributes_str or
        f'transcript_biotype "{bt}"' in attributes_str
        for bt in biotypes
    )

def detect_biotype_attribute(gtf_file: str) -> Optional[str]:
    """
    Detect which biotype attribute the GTF uses.
    Returns the attribute name or None if not found.
    """
    candidates = ['gene_type', 'gene_biotype', 'transcript_type', 'transcript_biotype']
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            for candidate in candidates:
                if candidate in line:
                    return candidate
    return None


def validate_gtf(file: str) -> None:
    """
    Validate that the file looks like a GTF:
    checks that non-comment lines have exactly 9 tab-separated fields.
    Raises ValueError if the file does not pass validation.
    """
    with open(file, "r") as f:
        for i, line in enumerate(f):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GTF: line {i+1} has {len(fields)} fields (expected 9).\\n"
                    f"  Line: {line.strip()[:120]}"
                )
            # Only check the first 100 data lines for speed
            if i > 100:
                break


def format_yaml_like(data, indent: int = 0) -> str:
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def filter_gtf_lncrna(gtf_in: str, filtered_gtf_out: str, biotypes: Set[str]) -> None:
    """
    Filter GTF file to extract only lncRNA entries.

    First pass:  identify all lncRNA transcript IDs from 'gene' and 'transcript' features.
    Second pass: write lines whose transcript_id is in the lncRNA set.
    """

    validate_gtf(gtf_in)

    # Warn if no biotype attribute is detected in the GTF
    biotype_attr = detect_biotype_attribute(gtf_in)
    if biotype_attr is None:
        logger.warning(
            "No biotype attribute found in GTF. Transcripts may not be correctly classified. "
            "Expected one of: gene_type, gene_biotype, transcript_type, transcript_biotype"
        )

    # --- First pass: collect lncRNA transcript IDs ---
    lncrna_transcript_ids: Set[str] = set()
    gene_count       = 0
    transcript_count = 0

    logger.info(f"First pass: identifying lncRNA entries in {gtf_in}")

    with open(gtf_in) as gtf:
        for line in gtf:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes   = fields[8]

            if not is_lncrna_biotype(attributes, biotypes):
                continue

            if feature_type == 'gene':
                # Count lncRNA genes (informational only)
                gene_count += 1

            elif feature_type == 'transcript':
                transcript_id = extract_attribute(attributes, 'transcript_id')
                if transcript_id:
                    lncrna_transcript_ids.add(transcript_id)
                    transcript_count += 1

    logger.info(f"Identified {gene_count} lncRNA genes")
    logger.info(f"Identified {transcript_count} lncRNA transcripts")

    if transcript_count == 0:
        logger.warning("No lncRNA entries found in GTF file!")

    # --- Second pass: write filtered GTF ---
    logger.info(f"Second pass: writing filtered GTF to {filtered_gtf_out}")

    # Features that belong to a transcript (resolved via transcript_id)
    TRANSCRIPT_FEATURES = {'exon', 'UTR', 'CCDS'}

    written_lines       = 0
    written_genes       = 0
    written_transcripts = 0
    written_features    = 0

    with open(gtf_in) as infile, open(filtered_gtf_out, "w") as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            if not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes   = fields[8]
            keep_line    = False

            if feature_type == 'gene':
                # Keep gene lines that are lncRNA biotype
                if is_lncrna_biotype(attributes, biotypes):
                    keep_line = True
                    written_genes += 1

            elif feature_type == 'transcript':
                # Keep transcript lines whose ID is in our lncRNA set
                transcript_id = extract_attribute(attributes, 'transcript_id')
                if transcript_id and transcript_id in lncrna_transcript_ids:
                    keep_line = True
                    written_transcripts += 1

            elif feature_type in TRANSCRIPT_FEATURES:
                # Keep sub-transcript features that belong to an lncRNA transcript
                transcript_id = extract_attribute(attributes, 'transcript_id')
                if transcript_id and transcript_id in lncrna_transcript_ids:
                    keep_line = True
                    written_features += 1

            # Note: CDS, start_codon, stop_codon are NOT included
            # because lncRNAs are non-coding by definition

            if keep_line:
                outfile.write(line)
                written_lines += 1

    if written_lines == 0:
        logger.warning("No lncRNA lines written to output file!")
    else:
        logger.info(f"Wrote {written_lines} total lines to {filtered_gtf_out}")
        logger.info(f"  - {written_genes} gene entries")
        logger.info(f"  - {written_transcripts} transcript entries")
        logger.info(f"  - {written_features} feature entries (exon, UTR, etc.)")


def main():
    parser = argparse.ArgumentParser(
        description='Filter GTF file to extract only lncRNA entries'
    )
    parser.add_argument('--gtf',    required=True, help='Input GTF file')
    parser.add_argument('--prefix', required=True, help='Output prefix')
    parser.add_argument('--biotypes', required=True, help='Comma-separated list of lncRNA biotypes (from params.lncrna_biotypes)')
    parser.add_argument('--log',    default=None,  help='Log file path (optional)')
    args = parser.parse_args()

    # Parse biotypes from comma-separated string into a set
    biotypes = set(bt.strip() for bt in args.biotypes.split(',') if bt.strip())
    if not biotypes:
        raise ValueError("--biotypes cannot be empty")

    # Set up file logging if requested
    log_file = args.log or f"{args.prefix}.filter.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(name)s - %(asctime)s %(levelname)s: %(message)s")
    )
    logger.addHandler(file_handler)

    output_file = f"{args.prefix}.lncrna.gtf"

    try:
        logger.info("Starting lncRNA GTF filtering")
        logger.info(f"Input GTF:  {args.gtf}")
        logger.info(f"Output GTF: {output_file}")
        logger.info(f"Biotypes:   {sorted(biotypes)}")

        filter_gtf_lncrna(args.gtf, output_file, biotypes)

        logger.info(f"Successfully created: {output_file}")
        print(f"Successfully created: {output_file}")

    except Exception as e:
        logger.error(f"Error filtering GTF: {e}")
        print(f"Error filtering GTF: {e}")
        raise


if __name__ == '__main__':
    import sys

    # Nextflow template substitution mode
    if len(sys.argv) == 1:
        sys.argv = [
            "filter_gtf_lncrna.py",
            "--gtf",      "$gtf",
            "--prefix",   "$gtf.baseName",
            "--biotypes", "$lncrna_biotypes",
        ]

    main()

    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
