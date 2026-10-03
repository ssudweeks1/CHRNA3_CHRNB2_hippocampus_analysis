# CHRNA3/CHRNB2 Hippocampal Transcript Co-detection Analysis

This repository contains the R analysis scripts used for the publicly available single-cell and single-nucleus RNA-sequencing analyses reported in:

**CHRNA3 and CHRNB2 transcript expression and co-detection in rat, mouse, and human hippocampus**

The analyses examine transcript detection and same-cell/same-nucleus co-detection of the nicotinic acetylcholine receptor subunit genes **Chrna3/Chrnb2** in hippocampal inhibitory neuronal populations.

## Analyses included

### Adult mouse

**Script:** `01_Allen_mouse_CHRNA3_CHRNB2_figure2_revision.R`

Analysis of the Allen Institute for Brain Science Mouse Whole Cortex and Hippocampus – 10x dataset.

The analysis:
- restricts cells to the hippocampal formation and GABAergic neuronal class;
- classifies cells as Chrna3-only, Chrnb2-only, double-positive, or negative for both transcripts;
- evaluates double-positive frequency across inhibitory transcriptomic classes;
- performs focused analyses of selected transcriptomic populations; and
- generates the analyses and figure corresponding to Figure 2.

### Mouse aging

**Script:** `Allen_Mouse_Aging_CHRNA3_CHRNB2_analysis.R`

Analysis of the Allen Institute for Brain Science Cellular and Molecular Characterization of the Aged Mouse Brain dataset.

The primary aging analysis:
- restricts cells to hippocampal GABAergic populations;
- uses the unbiased sampling arm as the prespecified primary cohort;
- classifies individual cells using a raw UMI count > 0 for transcript detection;
- treats the individual mouse as the biological replicate for age comparisons;
- compares young-adult and aged mice using exact donor-level permutation tests; and
- includes all-sampling and genotype-matched sensitivity analyses.

### Human hippocampus

**Script:** `Thompson_Human_CHRNA3_CHRNB2_Figure4.R`

Analysis of the anterior human hippocampal single-nucleus RNA-sequencing dataset reported by Thompson et al. (2025) and deposited under GEO accession **GSE264624**.

The analysis:
- restricts nuclei to inhibitory/GABAergic populations;
- uses unnormalized gene counts with a detection threshold of > 0;
- classifies nuclei as CHRNA3-only, CHRNB2-only, double-positive, or negative for both genes;
- evaluates double-positive frequency across broad inhibitory classes and finer transcriptomic populations;
- performs donor-aware analyses of the GABA.MGE CORT population; and
- generates the analyses and figures corresponding to Figures 4 and 5.

## Data sources

The underlying sequencing datasets are publicly available from the original data repositories and are not redistributed in this repository.

**Adult mouse:** Allen Institute for Brain Science, Mouse Whole Cortex and Hippocampus – 10x dataset (Yao et al., 2021).

**Mouse aging:** Allen Institute for Brain Science, Cellular and Molecular Characterization of the Aged Mouse Brain (10x-scRNAseq-aged-adult) dataset. The analysis used:
- Allen release manifest: `releases/20260711/manifest.json`
- metadata version: `Zeng-Aging-Mouse-10Xv3`, `20250131`
- raw expression matrix version: `20241130`
- taxonomy version: `Zeng-Aging-Mouse-WMB-taxonomy`, `20241130`

**Human:** Thompson et al. (2025), anterior human hippocampus single-nucleus transcriptomic dataset, GEO **GSE264624**.

## Transcript detection

For the single-cell and single-nucleus analyses, transcript/gene detection was defined as a raw count greater than zero.

## Reproducibility and biological replication

For the mouse aging analysis, individual mice were treated as biological replicates. Individual cells were used to calculate donor-level proportions but were not treated as independent biological replicates for age inference.

The human donor-aware analyses likewise preserve the individual donor as the biological replicate.

## Software

The publicly available dataset analyses were performed in:

- R 4.6.1
- RStudio 2026.08.0 Build 187
- ggplot2
- patchwork
- ggtext

## Citation

Please cite the manuscript and the original data sources when using these analyses.

## Contact

For questions regarding the analysis code, please contact the corresponding author.
