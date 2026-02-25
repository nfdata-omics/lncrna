#!/usr/bin/env python3
import re
import os
import platform


def has_pc(attrs: str) -> bool:
    return (
        'gene_biotype "protein_coding"' in attrs
        or 'gene_type "protein_coding"' in attrs
        or 'transcript_biotype "protein_coding"' in attrs
        or 'transcript_type "protein_coding"' in attrs
    )


def main():
    gtf_path = "$gtf"
    prefix = os.path.splitext(os.path.basename(gtf_path))[0]

    pc_genes = set()
    pc_txs = set()

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            attrs = parts[8]
            if feature == "gene" and has_pc(attrs):
                m = re.search(r'gene_id "([^"]+)"', attrs)
                if m:
                    pc_genes.add(m.group(1))
            elif feature == "transcript" and has_pc(attrs):
                m = re.search(r'transcript_id "([^"]+)"', attrs)
                if m:
                    pc_txs.add(m.group(1))

    out_path = f"{prefix}.protein_coding.gtf"
    with open(out_path, "w") as out, open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            gid_m = re.search(r'gene_id "([^"]+)"', attrs)
            tid_m = re.search(r'transcript_id "([^"]+)"', attrs)
            gid = gid_m.group(1) if gid_m else None
            tid = tid_m.group(1) if tid_m else None
            if (gid and gid in pc_genes) or (tid and tid in pc_txs):
                out.write(line)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {platform.python_version()}\\n')


if __name__ == "__main__":
    main()
