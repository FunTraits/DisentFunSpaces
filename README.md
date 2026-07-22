# Disentangling Functional Spaces

[![Data DOI](https://img.shields.io/badge/Data-10.5281%2Fzenodo.21487237-blue)](https://doi.org/10.5281/zenodo.21487237)
<!-- Code archive badge — replace XXXXXXX with the DOI Zenodo mints on your first GitHub release:
[![Code DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

This repository contains the code, figures, and result tables supporting the article:

**Title:** *Disentangling functional spaces: Toward a pluralistic view of ecological roles in conservation*
**Journal:** *Ecography*
**Authors:** Aurèle Toussaint, Pablo Tedesco, Gaël Grenouillet, Liis Kasari-Toussaint, Sébastien Brosse
**Affiliation:** CNRS, France

---

## 🗂️ Repository structure

```text
DisentFunSpaces/
├── scripts/                  # Analysis pipeline (sourced in numeric order by run_all.R)
│   ├── 00_START_GeneralScript.R          # Environment init, global options, shared utils
│   ├── 01_DATA_load_and_clean.R          # Preprocessing: traits, IUCN, imputation, PCoA, TPD
│   ├── 03_FIGURE_FUn.R                    # Main figures
│   ├── 03_FIGURE_FamilyFUn.R
│   ├── 03_FIGURE_IUCN.R
│   ├── 03_FIGURE_R2_aggregated_vs_disaggregated.R
│   ├── 03_FIGURE_spaceall.R
│   ├── 03_FIGURE_traitspaces.R
│   ├── 04_SUPP_*.R                        # Supplementary analyses (Beak, CFA, imputation, PCA robustness, …)
│   ├── 05_STAT_permutation_unique_species.R  # Permutation tests
│   ├── setup_environment.R               # Helper scripts (not run by run_all.R)
│
├── utils/                    # Custom functions and plotting helpers
│   ├── functions.R
│   ├── library.R
│   ├── generate_FRic_maps.R
│   ├── plot_MLD_map_cat.R
│   ├── plot_SES_map_with_latprofile.R
│   └── plot_diff_with_legend.R
│
├── results/                  # Outputs of the analyses
│   ├── figures/              # Main and supplementary figures (+ Framework.pptx)
│   └── tables/               # Summary and result tables (.csv, .rds)
│
├── data/                     # ⚠️ Not tracked (see Data availability) — raw/, processed/
│
├── run_all.R                 # Runs the full pipeline (sources numbered scripts sequentially)
├── DisentFunSpaces.Rproj     # RStudio project file
├── renv.lock.txt             # renv lockfile for environment reproducibility
├── LICENSE.txt               # License information
├── README.md                 # This file
└── .gitignore
```

> **Note:** `data/` are listed in `.gitignore` and are **not** included in the Git repository.

---

## 📦 Data availability

The `data/` folder (raw datasets, processed objects, and `metadata.xlsx`) is archived separately on **Zenodo** and is required to reproduce the analyses.

1. Download the data archive from Zenodo: [10.5281/zenodo.21487237](https://doi.org/10.5281/zenodo.21487237)
2. Unzip it at the root of the project so that the structure becomes:

```text
DisentFunSpaces/
├── data/
│   ├── raw/            # Original datasets (AVONET traits, BOTW ranges, phylogeny, IUCN, …)
│   ├── processed/      # Cleaned and formatted objects produced by 01_DATA_load_and_clean.R
```

---

## 🧪 Reproducing the analysis

This project uses [`renv`](https://rstudio.github.io/renv/) for dependency management.

1. Clone the repository:

   ```bash
   git clone https://github.com/FunTraits/DisentFunSpaces.git
   cd DisentFunSpaces
   ```

2. Add the `data/` folder downloaded from Zenodo (see **Data availability** above).

3. Open `DisentFunSpaces.Rproj` in RStudio, then restore the environment and run the full pipeline:

   ```r
   # renv.lock is provided as renv.lock.txt — rename it if needed:
   # file.rename("renv.lock.txt", "renv.lock")
   renv::restore()
   source("run_all.R")
   ```

   `run_all.R` sources every script in `scripts/` whose name matches `NN_*.R` (i.e. `00_`, `01_`, `03_`, `04_`, `05_`) in numeric order. Outputs are written to `results/figures/` and `results/tables/`.

---

## 📌 How to cite

If you use this code or data, please cite the article and the archives:

> Toussaint, A., Tedesco, P., Grenouillet, G., Kasari-Toussaint, L., & Brosse, S. Disentangling functional spaces: Toward a pluralistic view of ecological roles in conservation. *Ecography*. https://doi.org/<ARTICLE-DOI>

- **Data:** [10.5281/zenodo.21487237](https://doi.org/10.5281/zenodo.21487237)
- **Code:** *DOI minted by Zenodo on the first GitHub release — add it here.*

---

## 📄 License

See [`LICENSE.txt`](LICENSE.txt).

## ✒️ Contact

Aurèle Toussaint — CNRS, Toulouse, France — aurele.toussaint@cnrs.fr
