# Programs, origins and immunomodulatory functions of myeloid cells in glioma

[![Docker Build & Push](https://github.com/ilariag04/glioma-myeloid-programs-repro/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/ilariag04/glioma-myeloid-programs-repro/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Computational reproducibility workflow for the single-cell RNA-seq re-analysis of myeloid programs in glioma, adhering to **FAIR** (Findable, Accessible, Interoperable, Reusable) data principles.

---

## 📄 Study & Data Provenance

* **Study Title:** *Programs, origins and immunomodulatory functions of myeloid cells in glioma*
* **Accession:** Broad Institute Single Cell Portal (`SCP2389`)
* **Technology:** Seq-Well single-cell RNA-seq
* **Primary Target:** Microglial and macrophage functional phenotypes across WT and mutant brain tumors.

---

## 📁 Repository Architecture

```text
├── .github/
│   └── workflows/
│       └── docker-publish.yml    # Automated container deployment to GHCR
├── data/
│   └── download_SCP2389.sh       # Instructions and shell script for raw data retrieval
├── figures/                      # Output quality control, Elbow, and UMAP plots
├── scRNAseq/
│   └── 01_processing.r           # End-to-end QC, SCTransform v2, PCA, and Louvain clustering
├── Dockerfile                    # Container definition (R 4.x, Seurat v5, sctransform)
├── LICENSE                       # MIT License
└── README.md

```
Computational Methodology
1. Count Ingestion: Direct sparse matrix loading (matrix.mtx, barcodes.tsv, genes.tsv) into a Seurat object via Read10X.
2. Quality Control (QC): Elimination of damaged cells, empty droplets, and unresolved doublets:
   * Detected genes: `500 < nFeature_RNA < 6000`
   * UMI counts: `nCount_RNA > 1000`
   * Mitochondrial content: `percent.mt < 20%`
4. Variance Stabilization: Modeling expression with SCTransform (regularized negative binomial regression, vst.flavor = "v2"), regressing out mitochondrial contamination.
5. Dimensionality Reduction: Run on 50 principal components followed by Elbow Plot inspection.
6. Clustering: Louvain community detection (dims = 1:24, resolution = 0.3) projected into 2D via UMAP.
7. Memory Optimizations: Mitigation of Out-of-Memory (OOM) failures by clearing dense scale.data matrices and saving the final object uncompressed (saveRDS(..., compress = FALSE)).

Container EnvironmentThe analysis environment is containerized to guarantee identical library dependencies across systems.

Pull from GitHub Container Registry
docker pull ghcr.io/ilariag04/glioma-myeloid-programs-repro:latest
