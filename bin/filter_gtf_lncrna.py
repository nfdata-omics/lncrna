#!/usr/bin/env python3

"""
Filter GTF file to extract only lncRNA entries.
Designed to work after filter_gtf.py for validation.

Written by Karla Ruiz. Released under the MIT license.
"""

import logging
import argparse
import re
import statistics
from typing import Optional, Set

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("filter_gtf_lncrna")
logger.setLevel(logging.INFO)


def tab_delimited(file: str) -> float:
    """Check if file is tab-delimited and return median number of tabs."""
    with open(file, "r") as f:
        data = f.read(102400)
        tabs_per_line = [line.count("\t") for line in data.split("\n") if line.strip()]
        return statistics.median(tabs_per_line) if tabs_per_line else 0


def extract_attribute(attributes_str: str, attribute_name: str) -> Optional[str]:
    """
    Extract attribute value from GTF attributes string.

    Example: gene_id "ENSG00000000003"; gene_name "DDX11L1"
    Returns: ENSG00000000003
    """
    pattern = f'{attribute_name} "([^"]+)"'
    match = re.search(pattern, attributes_str)
    if match:
        return match.group(1)
    return None


def is_lncrna_biotype(attributes_str: str) -> bool:
    """Check if an attribute string contains lncRNA biotype."""
    return ('gene_type "lncRNA"' in attributes_str or
            'gene_biotype "lncRNA"' in attributes_str or
            'transcript_type "lncRNA"' in attributes_str or
            'transcript_biotype "lncRNA"' in attributes_str)


def filter_gtf_lncrna(gtf_in: str, filtered_gtf_out: str) -> None:
    """
    Filter GTF file to extract only lncRNA entries.

    First pass: Identify all lncRNA gene IDs and transcript IDs
    Second pass: Write only lines associated with lncRNA entries
    """

    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 9 tab-separated columns.")

    # Track lncRNA gene and transcript IDs
    lncrna_gene_ids: Set[str] = set()
    lncrna_transcript_ids: Set[str] = set()

    logger.info(f"First pass: identifying lncRNA entries in {gtf_in}")

    # First pass: collect lncRNA IDs
    gene_count = 0
    transcript_count = 0

    try:
        with open(gtf_in) as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                attributes = fields[8]

                # Process gene features
                if feature_type == 'gene':
                    if is_lncrna_biotype(attributes):
                        gene_id = extract_attribute(attributes, 'gene_id')
                        if gene_id:
                            lncrna_gene_ids.add(gene_id)
                            gene_count += 1

                # Process transcript features
                elif feature_type == 'transcript':
                    if is_lncrna_biotype(attributes):
                        transcript_id = extract_attribute(attributes, 'transcript_id')
                        gene_id = extract_attribute(attributes, 'gene_id')

                        if transcript_id:
                            lncrna_transcript_ids.add(transcript_id)
                            transcript_count += 1
                        if gene_id:
                            lncrna_gene_ids.add(gene_id)

        logger.info(f"Identified {gene_count} lncRNA genes")
        logger.info(f"Identified {transcript_count} lncRNA transcripts")

        if gene_count == 0 and transcript_count == 0:
            logger.warning("No lncRNA entries found in GTF file!")

        # Second pass: write filtered GTF
        logger.info(f"Second pass: writing filtered GTF to {filtered_gtf_out}")

        written_lines = 0
        written_genes = 0
        written_transcripts = 0
        written_features = 0

        with open(gtf_in) as infile, open(filtered_gtf_out, "w") as outfile:
            for line in infile:
                # Always keep header lines
                if line.startswith('#'):
                    outfile.write(line)
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                attributes = fields[8]

                keep_line = False

                # Keep lncRNA gene entries
                if feature_type == 'gene':
                    if is_lncrna_biotype(attributes):
                        keep_line = True
                        written_genes += 1

                # Keep lncRNA transcript entries
                elif feature_type == 'transcript':
                    if is_lncrna_biotype(attributes):
                        keep_line = True
                        written_transcripts += 1

                # Keep feature entries (exon, CDS, etc.) associated with lncRNA transcripts
                elif feature_type in ['exon', 'CDS', 'start_codon', 'stop_codon',
                                      'Selenocysteine', 'UTR', 'CCDS']:
                    transcript_id = extract_attribute(attributes, 'transcript_id')
                    if transcript_id and transcript_id in lncrna_transcript_ids:
                        keep_line = True
                        written_features += 1

                if keep_line:
                    outfile.write(line)
                    written_lines += 1

        if written_lines == 0:
            logger.warning("No lncRNA lines written to output file!")
        else:
            logger.info(f"Wrote {written_lines} total lines to {filtered_gtf_out}")
            logger.info(f"  - {written_genes} gene entries")
            logger.info(f"  - {written_transcripts} transcript entries")
            logger.info(f"  - {written_features} feature entries (exon, CDS, etc.)")

    except IOError as e:
        logger.error(f"File operation failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Filter GTF file to extract only lncRNA entries. "
                    "Designed to work after fasta_gtf_filter.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - filter annotated lncRNAs
  filter_gtf_lncrna.py --gtf genome.filtered.gtf --prefix output

  # With logging to file
  filter_gtf_lncrna.py --gtf genome.filtered.gtf --prefix lncrna \\
    --log_file filter.log
        """
    )

    parser.add_argument("--gtf", type=str, required=True,
                        help="Input GTF file (should be pre-filtered by fasta_gtf_filter.py)")
    parser.add_argument("--prefix", type=str, required=True,
                        help="Output prefix (will create PREFIX.lncrna.gtf)")
    parser.add_argument("--log_file", type=str, default=None,
                        help="Optional log file for detailed output")

    args = parser.parse_args()

    # Setup file logging if requested
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(name)s - %(asctime)s %(levelname)s: %(message)s")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    output_file = f"{args.prefix}.lncrna.gtf"

    try:
        logger.info(f"Starting lncRNA GTF filtering")
        logger.info(f"Input GTF: {args.gtf}")
        logger.info(f"Output GTF: {output_file}")

        filter_gtf_lncrna(args.gtf, output_file)

        logger.info(f"✓ Successfully created: {output_file}")
        print(f"✓ Successfully created: {output_file}")

    except Exception as e:
        logger.error(f"✗ Error filtering GTF: {e}")
        print(f"✗ Error filtering GTF: {e}")
        exit(1)


if __name__ == '__main__':
    main()
