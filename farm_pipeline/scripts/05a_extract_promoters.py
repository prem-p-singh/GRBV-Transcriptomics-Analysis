#!/usr/bin/env python
"""05a. Extract gene promoters (N bp upstream of TSS) from genome FASTA + GFF.

Outputs a FASTA with one record per gene, keyed by NCBI gene symbol (LOC…).
"""
import argparse, re, sys
from collections import defaultdict
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--genome", required=True)
ap.add_argument("--gff",    required=True)
ap.add_argument("--upstream_bp", type=int, default=1500)
ap.add_argument("--out", required=True)
args = ap.parse_args()

# Parse GFF to get gene records: chrom, strand, start, end, gene_name
print("[05a] loading genome…", file=sys.stderr)
genome = {r.id: r for r in SeqIO.parse(args.genome, "fasta")}

genes = {}   # gene_name -> (chrom, strand, start, end)
with open(args.gff) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        attr = p[8]
        m = re.search(r'gene=([^;]+)', attr) or re.search(r'Name=([^;]+)', attr)
        if not m:
            continue
        name = m.group(1)
        chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
        if chrom not in genome:
            continue
        genes[name] = (chrom, strand, start, end)

print(f"[05a] parsed {len(genes)} genes", file=sys.stderr)

with open(args.out, "w") as out:
    n_written = 0
    for name, (chrom, strand, start, end) in genes.items():
        seq_rec = genome[chrom]
        L = len(seq_rec.seq)
        if strand == "+":
            s = max(0, start - 1 - args.upstream_bp); e = max(0, start - 1)
            prom = seq_rec.seq[s:e]
        else:
            s = end; e = min(L, end + args.upstream_bp)
            prom = seq_rec.seq[s:e].reverse_complement()
        if len(prom) < 100:
            continue
        out.write(f">{name}\n{str(prom)}\n")
        n_written += 1
print(f"[05a] wrote {n_written} promoter sequences", file=sys.stderr)
