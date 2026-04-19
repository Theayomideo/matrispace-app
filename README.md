# MatriSpace: Identification and visualization of spatially-resolved ECM gene expression patterns in health and disease

[Ayomide Oshinjo](https://orcid.org/0000-0003-4303-0181) 1, [Daiqing Chen](https://orcid.org/0009-0003-5996-5181) 2, [Petar B. Petrov](https://orcid.org/0000-0001-5551-8032) 1,3, [Valerio Izzi](https://orcid.org/0000-0002-9960-4917) 1,3\*, [Alexandra Naba](https://orcid.org/0000-0002-4796-5614) 2,4\*

1. Faculty of Biochemistry and Molecular Medicine & Faculty of Medicine, BioIM Unit, University of Oulu, Oulu, FI-90014, Finland
2. Department of Physiology and Biophysics, University of Illinois Chicago, Chicago, IL 60612, USA
3. Infotech Institute, University of Oulu, Oulu, FI-90014, Finland
4. University of Illinois Cancer Center, Chicago, IL 60612, USA

\* equally-contributing corresponding authors

---

* MatriSpace features an intuitive graphical user interface implemented in R Shiny.
* The online version of MatriSpace, supporting both user uploads and a collection of 198 pre-processed datasets, is available and **ready-to-use** at http://matrinet.shinyapps.io/matrispace
* This repository contains the **offline version** of MatriSpace for local use with user-uploaded data (no dataset size restrictions). See [Installation](#installation) instructions below.
* The [MatriSpace-analyses](https://github.com/Izzilab/MatriSpace-analyses) repository includes examples of analyses performed using MatriSpace.
* If you use MatriSpace in your publications, please cite our preprint: **doi** [xxxxx](https://doi.org/xxxxx)

[![Badge](https://img.shields.io/badge/MatriSpace-Online-blue)](http://matrinet.shinyapps.io/matrispace)
[![Badge](https://img.shields.io/badge/Installation-info-green)](#installation)
[![Badge](https://img.shields.io/badge/Analysis-examples-orange)](https://github.com/Izzilab/MatriSpace-analyses)
[![Badge](https://img.shields.io/badge/Manuscript-bioRxiv-red)](https://doi.org/xxxxx)
[![Badge](https://img.shields.io/badge/Release-v1.0-green)](https://github.com/Theayomideo/matrispace-app/releases/tag/v1.0)

* Authors and maintainers: [Izzi Lab](https://www.oulu.fi/en/research-groups/izzi-group) and [Naba Lab for ECM Research](https://naba.lab.uic.edu/) (✉️ <matrisomeproject@gmail.com>)
* This work was supported by the following grants (green: Naba Lab; blue: Izzi Lab):

[![Badge](https://img.shields.io/badge/HuBMAP-U01HG012680-lightgreen)](https://commonfund.nih.gov/HuBMAP)
[![Badge](https://img.shields.io/badge/IMAT-R21CA261642-lightgreen)](https://www.cancer.gov/about-nci/organization/cssi/research/imat)
[![Badge](https://img.shields.io/badge/EU-CARES-lightblue)](https://www.cares-eu.org/)
[![Badge](https://img.shields.io/badge/CFF-2023--2024-lightblue)](https://syopasaatio.fi/)
[![Badge](https://img.shields.io/badge/DigiHealth-Oulu_University/Infotech-lightblue)](https://www.oulu.fi/en/research/creating-better-health-our-digital-health-knowhow)

## Motivation

The extracellular matrix (ECM) is a complex meshwork of proteins that serves as the master structural organizer of all tissues. Although spatial transcriptomics enables the study of ECM gene expression in its native tissue context, existing computational tools are designed primarily for intracellular expression variability and cell-cell signaling, and thus poorly capture the spatial architecture, regional specialization, and niche-specific deployment of ECM components.

MatriSpace is a computational framework designed to identify, quantify, and interpret spatial expression patterns of matrisome genes, ECM gene sets, and ECM-defined niches from spatial transcriptomics data. It focuses on the spatial organization, co-localization, and regional coordination of matrisome components, as well as their relationships with non-matrisome genes across entire tissue sections.

## Graphical Abstract

[placeholder]

## Workflow

MatriSpace operates through three analytical stages: **Data Input**, **Matrisome Profiling**, and **Feature Analysis**.

### Data Input

MatriSpace accepts Seurat objects (`.rds`) or SpatialExperiment objects (`.rds`), ensuring compatibility with both Seurat and Bioconductor ecosystems. Upon loading, uploaded datasets undergo automated preprocessing including gene symbol standardization, matrisome gene set scoring, and ECM niche classification.

The [online version](http://matrinet.shinyapps.io/matrispace) supports both user uploads (up to 1 GB) and a curated collection of 198 pre-processed 10x Visium datasets from public repositories (10x Genomics, HTAN, GEO, Zenodo), spanning 12 cancer types (n = 180 samples) and 10 healthy organ systems (n = 18 samples).

### Matrisome Profiling

All analyses are stratified by a user-selected annotation variable (e.g., cell type labels, clusters, or pathologist-defined regions). The matrisome profiling workflow computes expression scores for **6 matrisome categories** (Glycoproteins, Collagens, Proteoglycans, ECM Regulators, Secreted Factors, ECM-affiliated), **4 subcategories** (Perivascular, Hemostasis, Elastic fibers, Growth-factor binding), and **4 matrisome gene families** (Laminins, Matricellular proteins, Syndecans, Glypicans).

Results are presented as spatial distribution maps and hotspot maps, each accompanied by spatial autocorrelation statistics.

**ECM niche classification.** Each spot is classified into Interstitial or Basement membrane ECM niches. Niche assignments and per-spot niche scores are visualized on an interactive spatial viewer.

**Spatial ligand-receptor co-expression.** ECM-focused interaction pairs from [MatriComDB](https://github.com/Izzilab/MatriCom) are evaluated by computing spatial co-expression scores within and across ECM niches, and visualized as interactive heatmaps and volcano plots.

### Feature Analysis

Users may interrogate specific features using a two-level selection system. The primary feature can be a matrisome gene or gene set, while the optional secondary feature may be any detected gene or matrisome gene set. Outputs include spatial feature maps, expression distributions across annotation groups, co-expression blend plots, and local indicators of spatial association (LISA) maps.

All results are exportable as plots and tables, and processed Seurat objects can be downloaded for further analysis.

## Installation

MatriSpace is distributed as a Docker image and an R package. The offline version supports user-uploaded data only; the pre-processed dataset collection is available exclusively through the [online version](http://matrinet.shinyapps.io/matrispace).

### Docker (recommended)

Requires [Docker Desktop](https://docs.docker.com/get-started/get-docker/) to be installed and running.

```bash
docker pull ghcr.io/theayomideo/matrispace-app:main
docker run -p 3838:3838 ghcr.io/theayomideo/matrispace-app:main
```

Once the container is running, open http://localhost:3838 in your browser.

### R Package

Requires [R](https://cran.r-project.org/) (>= 4.3.3).

```r
install.packages("pak")
pak::pkg_install("Theayomideo/matrispace-app")
matrispace.app::run_app()
```

## See Also

* [The Matrisome Project](https://sites.google.com/uic.edu/matrisome/home): an open-access resource sharing protocols, tools, and datasets to support ECM research.
* [MatriCom](https://github.com/Izzilab/MatriCom) and [MatriComDB](https://github.com/Izzilab/MatriCom): a tool and curated database for inferring ECM communication systems from scRNA-seq data.
* [Naba Lab for ECM Research](https://naba.lab.uic.edu/) ([@GitHub](https://github.com/Matrisome/)) and [Izzi Lab](https://www.oulu.fi/en/research-groups/izzi-group) ([@GitHub](https://github.com/izzilab))
