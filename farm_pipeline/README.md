# GRBV Extended Analysis Pipeline — FARM HPC

Snakemake pipeline running 10 analyses on top of your existing GRBV transcriptomics data, designed for UC Davis FARM (`low` partition, 20 GB quota, `publicgrp` account).

## Analyses included

| # | Rule | Tool | What it adds |
|---|------|------|--------------|
| 1 | `grn_inference` | GENIE3 | Directed TF → target edges (regulators = PlantTFDB Vvi) |
| 2 | `module_preservation` | WGCNA::modulePreservation | Tests whether RB− co-expression modules break down in RB+ |
| 3 | `subtype_discovery` | Consensus NMF | Transcriptional subtypes of infection in RB+ samples |
| 4 | `diffcoex` | DGCA | Gene pairs with significantly changed correlation RB+ vs RB− |
| 5 | `motif_enrichment` | FIMO + JASPAR plants | Cis-motifs enriched in each WGCNA module's promoters |
| 6 | `pseudotime` | slingshot + tradeSeq | Infection-age trajectory and temporally variable genes |
| 7 | `mofa_integration` | MOFA+ | Latent factors across expression + titer + weather |
| 8 | `hima_mediation` | HIMA | Multi-mediator Temperature → Genes → Titer |
| 9 | `defense_repertoire` | keyword + DE | Curated NLR/PRR/PR/SA/JA/ET panel |
| 10 | `subtype_characterize` | t-test + χ² | DE + clinical associations for subtypes from rule 3 |

## One-time setup on FARM

From your Mac:

```bash
cd ~/Data_Analysis/farm_pipeline
bash push.sh                               # rsync pipeline + inputs to farm:~/grbv_extended/
ssh farm "cd ~/grbv_extended && bash setup_farm.sh"   # conda env + genome + TF + motif downloads (~2–3 GB)
```

`setup_farm.sh` creates the `grbv-extended` conda env and downloads:
- NCBI V. vinifera reference genome + GFF (`GCF_030704535.1`)
- PlantTFDB Vvi TF list
- JASPAR2024 CORE plants motif database

## Run

```bash
ssh farm "cd ~/grbv_extended && bash run_farm.sh --cores 40"
# or dry-run first:
ssh farm "cd ~/grbv_extended && bash run_farm.sh -n"
```

All rules submit as separate SLURM jobs via the `profiles/farm/` executor profile. They are preempt-safe (`--requeue=true`). Expected walltimes per rule:

| Rule | Walltime | CPUs | Mem |
|---|---|---|---|
| 01 GRN | 6 h | 16 | 32 G |
| 02 mod-pres | 3 h | 8 | 16 G |
| 03 NMF | 1.5 h | 4 | 8 G |
| 04 DGCA | 2 h | 8 | 16 G |
| 05 motif (05a+05b) | 1.5 h | 4 | 8 G |
| 06 pseudotime | 1 h | 4 | 8 G |
| 07 MOFA | 2 h | 4 | 16 G |
| 08 HIMA | 3 h | 8 | 16 G |
| 09 defense | 30 m | 2 | 4 G |
| 10 subtype char. | 30 m | 2 | 4 G |

## Pull results back

```bash
cd ~/Data_Analysis/farm_pipeline
bash pull.sh                               # rsync results/ + logs/ from FARM
```

## Layout

```
farm_pipeline/
├── Snakefile                    orchestrator (10 rules)
├── config/config.yaml           thresholds / URLs
├── profiles/farm/config.yaml    SLURM executor (low, publicgrp, --requeue)
├── envs/grbv-extended.yaml      conda env spec
├── scripts/                     10 analysis scripts (R + Python)
├── data/                        inputs staged by push.sh
├── resources/                   downloads populated by setup_farm.sh
├── results/                     outputs (rsync'd back via pull.sh)
├── logs/                        per-rule stderr/stdout
├── setup_farm.sh                one-time bootstrap on FARM
├── push.sh                      Mac → FARM sync
├── pull.sh                      FARM → Mac sync
└── run_farm.sh                  snakemake launcher
```

## Notes / gotchas

- **Disk (20 GB cap):** genome + GFF + motifs ≈ 2–3 GB. Data inputs ≈ 25 MB. Intermediate outputs are CSV/PNG only (no BAM). Total usage should stay <6 GB.
- **Preemption on `low`:** every rule sets `requeue=true` and the SLURM account/partition/QOS from `config.yaml`.
- **MEs.csv mismatch:** your existing `MEs.csv` uses different sample IDs (`16H1`…) than `counts_matrix.txt` (`S01`…). All analyses that need module eigengenes recompute them from counts — no dependency on the mismatched file.
- **Stage spelling:** viral-titer file uses `Pre_véraison` / `Véraison` (accented); scripts map these to the unaccented `Pre_veraison` used in `sample_metadata.txt`.
- **HIMA exposure is per-season temperature** (not per-sample) because daily weather files are not joined to each sample's Julian day. If you have an exact sampling date per sample, swap the `agg_tavg()` call for a window-aggregate around that date.
- **Module names in `Module_genes.xlsx`** must match sheet names; the motif enrichment script reads all sheets and assumes column 1 contains gene IDs.

## Rerunning a single rule

```bash
ssh farm "cd ~/grbv_extended && bash run_farm.sh -R module_preservation"
```

## Troubleshooting

| Issue | Fix |
|---|---|
| `conda: command not found` in SLURM job | rule wrapper does `source /etc/profile.d/modules.sh; module load conda/base` — if still failing, check `logs/<rule>.log` for earlier error |
| `fimo: not found` | `meme` not installed in env. Activate env and run `meme -version` to confirm |
| HIMA runs out of memory | lower `hima.top_var` in `config/config.yaml` |
| GENIE3 too slow | lower `grn.n_trees` (default 1000 → 500) |
| Job preempted mid-run | Snakemake + `--requeue=true` will auto-resubmit; `restart-times: 2` in the profile retries twice |
