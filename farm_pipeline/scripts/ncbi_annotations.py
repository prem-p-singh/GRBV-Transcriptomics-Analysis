#!/usr/bin/env python
"""Fix #11: pull NCBI gene annotations for every LOC identifier named in the
manuscript (headline biomarkers, mediators, hub TFs, motif regulators).

Uses the already-downloaded NCBI V. vinifera GFF (genome.gff) to avoid
hitting the NCBI E-utilities API. Extracts gene names, synonyms, and
product descriptions for each LOC ID.
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gff", default="/Users/prempratapsingh/Data_Analysis/farm_pipeline/resources/genome/genome.gff")
ap.add_argument("--out", default="/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/annotations/gene_annotations.csv")
args = ap.parse_args()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)

# Every LOC gene named in the manuscript, integrated-analysis, extended-report,
# sensitivity-report, or VIGS plan. Update as the paper evolves.
KEY_GENES = [
    # Host biomarkers (paper Fig 20-22)
    "LOC100854143", "LOC100256964", "LOC104881857", "LOC100268011", "LOC100257913",
    # Top subtype DE (paper Fig 10, biomarker_loso_cv)
    "LOC100260258", "LOC100257403", "LOC100261580", "LOC109121469", "LOC100258410",
    "LOC100241782", "LOC100257854", "LOC100252646", "LOC100254233", "LOC100262782",
    "LOC100255325", "LOC100241410",
    # DiffCoEx hubs (rule 04)
    "LOC100253651", "LOC100258046", "LOC100245905", "LOC100255991", "LOC100242338",
    "LOC100248652", "LOC100260313", "LOC100261077", "LOC100246293", "LOC100855168",
    "LOC100854991",
    # GRN hub TFs (rule 01)
    "LOC100241332", "LOC100255815", "LOC100251326", "LOC100262112", "LOC100250388",
    "LOC100255488", "LOC100263995", "LOC100248184", "LOC100248013", "LOC100243201",
    "LOC100264963", "LOC100232838", "LOC100267536", "LOC100265549",
    # Pseudotime top genes (rule 06)
    "LOC100260288", "LOC100262729", "LOC100245911", "LOC100263489", "LOC100264626",
    "LOC100267804", "LOC100264738", "LOC100244103", "LOC109122806", "LOC100242302",
    "LOC100232901", "LOC100261426",
    # Robust mediators (HIMA + main mediation sensitivity)
    "LOC100266309", "LOC100232974", "LOC100249249", "LOC132254135",
    "LOC100253072", "LOC100259589", "LOC100260566", "LOC100257676", "LOC100852723",
    "LOC100853778", "LOC100261367", "LOC100854026", "LOC100247528",
    # HIMA SIS pool (top few)
    "LOC100264214", "LOC100255569", "LOC100264704", "LOC132253503",
    # ML stability (biomarker_loso_cv)
    "LOC100261750", "LOC100241372", "LOC100258268", "LOC100249883", "LOC100233114",
    "LOC100244812", "LOC100267121",
    # Symbols & named genes for completeness
    "N761_gp1", "N761_gp2", "N761_gp3",
    "GRIP61", "GRIP22", "RD22-C", "CBF4", "PIP1;2",
]
KEY_GENES = sorted(set(KEY_GENES))
print(f"[ann] resolving {len(KEY_GENES)} gene identifiers", file=sys.stderr)

# Parse GFF for gene records
annotations = {}
with open(args.gff) as f:
    for ln in f:
        if ln.startswith("#") or not ln.strip(): continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        attr = p[8]
        name = None
        m = re.search(r"gene=([^;]+)", attr)
        if m: name = m.group(1)
        syn = []
        ms = re.search(r"gene_synonym=([^;]+)", attr)
        if ms: syn = [s.strip() for s in ms.group(1).split(",")]
        prod = None
        mp = re.search(r"Note=([^;]+)", attr) or re.search(r"description=([^;]+)", attr)
        if mp: prod = mp.group(1)
        biotype = None
        mb = re.search(r"gene_biotype=([^;]+)", attr)
        if mb: biotype = mb.group(1)
        gid = None
        mg = re.search(r"Dbxref=GeneID:(\d+)", attr)
        if mg: gid = mg.group(1)
        chrom = p[0]; start = int(p[3]); stop = int(p[4])
        rec = dict(name=name, synonyms=";".join(syn), product=prod,
                   biotype=biotype, entrez_id=gid,
                   seqid=chrom, start=start, stop=stop)
        if name: annotations[name] = rec
        for s in syn:
            annotations.setdefault(s, rec)

rows = []
for g in KEY_GENES:
    r = annotations.get(g, {})
    rows.append({
        "gene_id": g,
        "ncbi_gene_name": r.get("name", ""),
        "synonyms": r.get("synonyms", ""),
        "product": r.get("product", ""),
        "biotype": r.get("biotype", ""),
        "entrez_id": r.get("entrez_id", ""),
        "seqid": r.get("seqid", ""),
        "start": r.get("start", ""),
        "stop": r.get("stop", ""),
        "resolved": bool(r),
    })
df = pd.DataFrame(rows)
df.to_csv(args.out, index=False)
resolved = df["resolved"].sum()
print(f"[ann] resolved {resolved}/{len(df)} gene identifiers; wrote {args.out}", file=sys.stderr)
print("\nUnresolved (likely non-canonical names or viral genes):", file=sys.stderr)
for g in df[~df["resolved"]]["gene_id"]:
    print(f"   {g}", file=sys.stderr)
