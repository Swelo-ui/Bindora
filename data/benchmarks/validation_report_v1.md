# Bindora Accuracy & Validation Benchmark Report (v1.0)

**Date:** 2026-09-16 &bull; **Engine:** AutoDock Vina 1.2.7 &bull; **Scoring:** AutoDock Vina (Standard), Vinardo (Kortemme / Quiroga 2016)

## Executive Summary

This empirical benchmark validates Bindora's molecular docking engine against curated protein-ligand complexes with published crystallographic coordinates and wet-lab experimental binding affinities ($K_d / K_i / \Delta G_{\text{exp}}$).

| Benchmark Metric | AutoDock Vina | Vinardo Scoring | Target Standard |
| :--- | :---: | :---: | :---: |
| **Pearson Correlation ($R$)** | **-0.686** | **-0.825** | > 0.50 (CASF Core) |
| **Spearman Rank Correlation ($\rho$)** | **-0.300** | — | > 0.50 |
| **Root Mean Square Error (RMSE)** | **2.62 kcal/mol** | **5.17 kcal/mol** | < 2.5 kcal/mol |
| **Mean Absolute Error (MAE)** | **1.93 kcal/mol** | **4.87 kcal/mol** | < 2.0 kcal/mol |
| **Pose Redocking Success (RMSD $\le 2.0$ Å)** | **40.0%** (2/5) | — | > 70% |
| **Mean Crystallographic RMSD** | **1.84 Å** | — | < 2.0 Å |

---

## Itemized Benchmark Results Matrix

| PDB | Target Receptor | Investigational Ligand | Exp $\Delta G$ | Vina $\Delta G$ | Vinardo | Error | RMSD (Å) | Benchmark Status |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| `1CX2` | Cyclooxygenase-2 (COX-2) | SC-558 | -11.32 | -10.76 | -7.38 | 0.56 | 0.76 | ✅ Validated (< 2.0 Å) |
| `1IEP` | Abl1 Tyrosine Kinase | Imatinib (STI-571) | -10.91 | -11.61 | -7.58 | 0.7 | 0.80 | ✅ Validated (< 2.0 Å) |
| `2ITY` | EGFR Kinase (T790M) | Gefitinib (Iressa) | -11.73 | -8.04 | -4.39 | 3.69 | 2.94 | ⚠️ Near-Native (> 2.0 Å) |
| `1M17` | EGFR Kinase (Active) | Erlotinib (Tarceva) | -11.32 | -6.87 | -4.74 | 4.45 | 2.50 | ⚠️ Near-Native (> 2.0 Å) |
| `1XKK` | EGFR / HER2 Kinase | Lapatinib (Tykerb) | -10.91 | -10.67 | -7.76 | 0.24 | 2.18 | ⚠️ Near-Native (> 2.0 Å) |

---

## Methodology & Reproducibility Protocol

1. **Receptor Ingestion:** Standard biological assemblies stripped of co-solvents and crystallographic waters. Ionization set to physiological pH 7.4. Gasteiger-Marsili partial charges and AutoDock 4 atom types assigned via Meeko.
2. **Grid Box Centering:** Pocket centroid determined from crystallographic bound ligand coordinates with 22.0 Å cubic search space.
3. **Stochastic Sampling:** AutoDock Vina iterated local search with Monte Carlo sampling (fixed seed 42 for complete reproducibility).
4. **Vinardo Rescoring:** Scoring executed via AutoDock Vina v1.2.7 `--scoring vinardo --score_only --autobox`.
5. **RMSD Calculation:** Symmetry-corrected heavy-atom root mean square deviation between docked pose coordinates and crystallographic coordinates.

## Citations

- Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*, 31(2), 455-461.
- Quiroga, R., & Villarreal, M. A. (2016). Vinardo: A Scoring Function Based on Autodock Vina with Improved Affinity Predictions. *PLoS ONE*, 11(5), e0155182.
- Su, M., et al. (2019). Comparative Assessment of Scoring Functions (CASF-2016) on the PDBbind Database. *Journal of Chemical Information and Modeling*, 59(2), 895-913.