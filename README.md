<div align="center">
  <h1>Grapevine Red Blotch Disease Pathogenesis in Ripening Berries</h1>
  <p><strong>Stage-resolved viral <em>C2</em> activation, host pathway disruption, and ripening delay across 44 field-grown Merlot vines</strong></p>
  <p>
    <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-1F8A84" alt="MIT License" /></a>
    <img src="https://img.shields.io/badge/3'_TAG--seq-327_Libraries-9E2A47" alt="327 TAG-seq libraries" />
    <img src="https://img.shields.io/badge/Petiole_dPCR-267_Assays-1F8A84" alt="267 dPCR assays" />
    <img src="https://img.shields.io/badge/Phenology-MPhS_r%3D0.986-54657E" alt="MPhS validated" />
    <img src="https://img.shields.io/badge/Status-In_Preparation-B8862D" alt="Manuscript in preparation" />
  </p>
</div>

<p align="center">
  <img src="figures/readme/graphical_abstract.svg" alt="GRBV disease pathogenesis across grape berry ripening stages" width="100%" />
</p>

> [!NOTE]
> **Pre-publication repository:** Because the manuscript is currently in preparation, unpublished manuscript files, main publication figures, count matrices, and raw sequencing reads are withheld until publication. This repository provides the reproducible statistical analysis code and high-level study overview.

## Project at a glance

- **Early antiviral silencing without cell death:** Around veraison, infected berries induce *RDR1*, *DCL2*, *NAC62*, *NTM1-like 9*, and *RSL1* while keeping all `160` expressed berry NLR receptors restrained.
- **Post-veraison viral *C2* switch & proteostasis collapse:** After veraison, complementary-sense *C2* rises (`+1.92`, `+0.95` log2), covarying locally with `681` host genes as *DCL2* subsides and cell-wall (*XTH7*), ribosome, 26S proteasome, and *sHSP* chaperones (`30/34` genes) drop alongside *metacaspase-9* induction.
- **Delayed ripening vs. sugar loss:** Foliar viral DNA load (`267` petiole dPCR chips; *Singh et al., 2026*) predicts berry maturation delay (`P = 0.0049`), while early-season leaf viral burden precedes severe harvest sugar deficits (`−8.9 to −10.9 °Bx`).

---

## Disrupted host pathways & viral dynamics

<p align="center">
  <img src="figures/readme/host_pathways.svg" alt="Major host factors and pathways disrupted by GRBV beyond maturation delay" width="100%" />
</p>

<p align="center">
  <img src="figures/readme/viral_host_dynamics.svg" alt="Leaf-petiole viral DNA, maturation delay vs must Brix loss, and berry viral ORF dynamics" width="100%" />
</p>

---

## Field cohort & analytical evolution

<p align="center">
  <img src="figures/readme/study_overview.svg" alt="44 Merlot vines sampled across 4 ripening stages and 2 seasons" width="100%" />
</p>

<p align="center">
  <img src="figures/readme/analysis_evolution.svg" alt="Analytical evolution from Era 1 screening to Era 8 pathology synthesis" width="100%" />
</p>

---

## Repository & reproduction

```text
.
├── Current_Analysis/
│   ├── simplified_reanalysis_2026-11/   # Mixed-model (dream), zenith, permutation & petiole dPCR scripts
│   └── mphs_external_check_2026-11/     # External Molecular Phenology Scale (MPhS v0.1.4) validation
├── figures/readme/                      # High-level overview schematics for this README
├── CITATION.cff
└── LICENSE
```

```bash
Rscript Current_Analysis/simplified_reanalysis_2026-11/01_data_audit.R
Rscript Current_Analysis/simplified_reanalysis_2026-11/02_developmental_score.R
Rscript Current_Analysis/simplified_reanalysis_2026-11/04_adjusted_dream.R
Rscript Current_Analysis/simplified_reanalysis_2026-11/10_petiole_dpcr_integration.R
python3 Current_Analysis/mphs_external_check_2026-11/verify_mphs.py
```

---

## Citation & acknowledgments

Code is [MIT licensed](LICENSE) ([`CITATION.cff`](CITATION.cff)); raw 3' TAG-seq reads and manuscript files will be released upon publication. Companion petiole dPCR dataset: *Singh et al. (2026), Arch. Virol. 171:183*.

In memoriam **Dr. Anita Oberholster** (UC Davis Viticulture & Enology). Conducted with **Dr. Mysore R. Sudarshana** (USDA-ARS / UC Davis Plant Pathology), **Dr. David E. Block**, and **Bradley N. Jenner**, supported by CDFA, USDA, and UC Davis FARM HPC. Maintained by **Prem Pratap Singh** (UC Davis).
