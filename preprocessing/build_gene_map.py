"""Build canonical (gene_id, gene_name) TSV from an Ensembl GTF.

Reads gene-feature lines, extracts gene_id and gene_name attributes (stripping
ENSG version suffix if present), dedupes, and writes a TSV. Used by
get_scimilarity_embeddings.py and get_scvi_embeddings.py for ENSEMBL<->symbol
mapping without BioMart calls.

Usage:
    python preprocessing/build_gene_map.py \\
        --gtf raw_data/reference/Homo_sapiens.GRCh38.114.gtf.gz \\
        --out data/gene_map/gene_map_grch38_v114.tsv
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def parse_gene_lines(gtf_path):
    opener = gzip.open if gtf_path.endswith(".gz") else open
    n_gene_lines = 0
    n_with_name = 0
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            n_gene_lines += 1
            attrs = dict(ATTR_RE.findall(fields[8]))
            gene_id = attrs.get("gene_id", "").split(".")[0]
            gene_name = attrs.get("gene_name", "").strip()
            if not gene_id or not gene_name:
                continue
            n_with_name += 1
            yield gene_id, gene_name
    print(f"GTF gene lines: {n_gene_lines}, with both id+name: {n_with_name}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", required=True, help="Path to (gzipped) Ensembl GTF")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    id_to_names = defaultdict(set)
    name_counts = defaultdict(int)
    for gid, gname in parse_gene_lines(args.gtf):
        id_to_names[gid].add(gname)
        name_counts[gname] += 1

    n_id_conflicts = sum(1 for names in id_to_names.values() if len(names) > 1)
    n_name_dups = sum(1 for c in name_counts.values() if c > 1)
    print(f"gene_ids with >1 distinct gene_name: {n_id_conflicts}", file=sys.stderr)
    print(f"gene_names appearing under >1 gene_id: {n_name_dups}", file=sys.stderr)
    assert n_id_conflicts == 0, f"{n_id_conflicts} gene_ids map to >1 name — resolve before writing"

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    n_rows = 0
    with open(args.out, "w") as f:
        f.write("gene_id\tgene_name\n")
        for gid in sorted(id_to_names):
            for gname in sorted(id_to_names[gid]):
                f.write(f"{gid}\t{gname}\n")
                n_rows += 1
    print(f"Wrote {n_rows} rows ({len(id_to_names)} unique gene_ids) -> {args.out}")


if __name__ == "__main__":
    main()
