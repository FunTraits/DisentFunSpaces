# The hidden cost of a single index: dimension-explicit functional vulnerability across the global avifauna

This repository contains the data, code, and figures supporting the article:

---

## 🗂️ Repository structure

```textdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
DisentanglingFunctionalSpaces/
├── data/               # Raw and processed data
│   ├── raw/            # Original datasets (trait data, IUCN, etc.)
│   ├── processed/      # Cleaned and formatted data used in the analyses
│   └── metadata.csv    # Description and sources of all datasets
├── scripts/            # Scripts used to run the analyses
├── figures/            # Outputs of the analyses
├── tables/            # Outputs of the analyses
├── utils/              # Custom functions used in the analyses
├── README.md           # This file
├── LICENSE             # License information (default: CC-BY 4.0)
├── .gitignore          # Files to exclude from version control
└── renv.lock           # Optional: for environment reproducibility

data are available in this repo: 

## 🧪 Reproducing the analysis

This project uses [`renv`](https://rstudio.github.io/renv/) for dependency management.  
To reproduce the exact package versions used in the analysis:

1. Clone the repository:

```bash
git clone https://github.com/AureleToussaint/FUSE_Bird.git

