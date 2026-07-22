# Disentangling Functional Spaces

This repository contains the data, code, and figures supporting the article:

**Title:** *Disentangling functional spaces: Toward a pluralistic view of ecological roles in conservation*  
**Journal:** *Ecography* 
**Author:** Aurele Toussaint, Pablo Tedesco, Gael Grenouillet, Liis Kasari-Toussaint, Sébastien Brosse
**Affiliation:** CNRS, France

---

## 🗂️ Repository structure

```textdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
DisentanglingFunctionalSpaces/
├── data/               # Raw and processed data
│   ├── raw/            # Original datasets (trait data, IUCN, etc.)
│   ├── processed/      # Cleaned and formatted data used in the analyses
│   └── metadata.csv    # Description and sources of all datasets
│
├── scripts/            # Scripts used to run the analyses
│   ├── 00_GeneralScript.R
│   ├── 01_load_and_clean.R
│   ├── 02_analyses_traitspaces.R
│   ├── 03_analyses_clades.R
│   └── 04_analyses_biogeo.R
│
├── results/            # Outputs of the analyses
│   ├── figures/        # Main and supplementary figures
│   ├── tables/         # Summary tables and results
│
├── utils/              # Custom functions used in the analyses
│   └── functions.R
│   └── library.R
│
├── README.md           # This file
├── LICENSE             # License information 
├── .gitignore          # Files to exclude from version control
└── renv.lock           # For environment reproducibility


## 🧪 Reproducing the analysis

This project uses [`renv`](https://rstudio.github.io/renv/) for dependency management.  
To reproduce the exact package versions used in the analysis:

1. Clone the repository:

```bash
git clone https://github.com/AureleToussaint/DisentanglingFunctionalSpaces.git