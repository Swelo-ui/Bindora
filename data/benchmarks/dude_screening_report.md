# Bindora DUD-E Virtual Screening Report (v2.0)

**Date:** 2026-09-19 | Subset: single_vegfr2 | Exh: 4

> **Mandatory Caveat (always reported alongside DUD-E results):**
> DUD-E decoys are property-matched (MW, cLogP, HBA, HBD, rotatable bonds) but NOT topologically diversified from actives. Some decoys carry structural similarity to actives (analogue bias), inflating enrichment scores for some methods. This caveat applies to all published DUD-E benchmarks (Mysinger et al. 2012). Results must be interpreted alongside pose-accuracy data (CASF-2016 redocking), not in isolation.

## Aggregate Statistics

| Metric | Value |
| :--- | :---: |
| Targets attempted | 1 |
| Targets successful | 1 |
| Mean ROC-AUC | 0.84 |
| Mean EF1% | 20.0 |
| Mean EF5% | 4.0 |
| Mean EF10% | 2.0 |

## Per-Target Results (Full Distribution)

| Target | Name | Data Source | Actives | Decoys | ROC-AUC | EF1% | EF5% | EF10% | Status |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| vegfr2 | VEGFR2 | ChEMBL | 5 | 10 | 0.840 | 20.00 | 4.00 | 2.00 | Good (AUC >= 0.7) |

## Data Sources & References

- **DUD-E**: Mysinger MM et al. J. Med. Chem. 2012, 55(14), 6582. DOI: 10.1021/jm300687e
- **ChEMBL** (fallback when DUD-E server unavailable): Mendez D et al. Nucleic Acids Res 2019, 47(D1):D930-D940. DOI: 10.1093/nar/gky1075
  - Actives: binding IC50 <= 1 µM (most potent first)
  - Inactives (decoy proxy): binding IC50 >= 50 µM