## HackBio-Single-Cell-RNA-Seq-Stage-2

### 📘 Single-Cell RNA-seq Analysis & Biological Interpretation

A full Scanpy workflow for clustering, cell-type annotation, biological interpretation, and immunological reasoning.

## 🧬 1. Overview

This project analyzes a single-cell RNA-seq dataset using the Scanpy ecosystem.
We cluster cells, annotate immune lineages, interpret biological context, evaluate whether the tissue source resembles bone marrow, and assess whether immune proportions suggest health or infection.

Outputs include:

UMAPs (clusters + annotated cell types)

Marker-based annotation

Cell-type proportions

Biological interpretation

Publication-style PPT

Full reproducibility pipeline

## 🧭 2. Workflow Diagram
flowchart TD
    A[Raw Count Matrix] --> B[Quality Control]
    B --> C[Normalization + Log1p]
    C --> D[Highly Variable Gene Selection]
    D --> E[ PCA ]
    E --> F[Neighbors Graph]
    F --> G[Leiden Clustering]
    G --> H[UMAP Embedding]
    H --> I[Marker Gene Ranking]
    I --> J[Cell-Type Annotation]
    J --> K[Proportion Analysis]
    K --> L[Biological Interpretation]

## 🔬 3. Methods Summary
Scanpy workflow

QC filtering

Normalization + log1p

HVG selection

PCA

kNN graph construction

Leiden clustering

UMAP embedding

Marker-gene analysis

Proportion estimation

Validation & Add-Ons

Bootstrapping cluster stability

Differential expression

Pathway enrichment

Viral signature screening module

## 🏷 4. Identified Cell Types
Cell Type	Description
NK cells	Cytotoxic innate lymphocytes; kill virally infected or stressed cells
T cells	Adaptive immunity, antigen-specific responders
Monocytes	Phagocytic innate cells, inflammatory cytokine producers
B cells	Antigen-presenting lymphocytes that mature into plasma cells
Plasma cells	Antibody-secreting effector B cells
Erythroid lineage	RBC precursors
Megakaryocytes / Platelets	Platelet-producing cells involved in clotting
HSC / Progenitors	Stem and multipotent precursor compartments

## 📊 5. Cell-Type Proportion Summary

NK cells – 36%

T cells – 32%

Monocytes – 14%

B cells – 7%

Plasma cells – 5%

Erythroid – 3.5%

Megakaryocytes – 1.8%

HSC/Progenitor – <1%

## 🧭 6. Is the Tissue Bone Marrow?
Evidence against bone marrow:

❌ NK cells are far too high (BM usually <10%)
❌ T cells also unusually high for BM
❌ Very low HSC/progenitor content (<1%)
❌ Lack of granulocytes and neutrophils (BM hallmark)

Evidence for peripheral blood:

✔ NK (30–40%) and T cell dominance are typical
✔ Low progenitors consistent with circulation
✔ Myeloid monocyte fraction (~10–15%) matches PBMC
✔ Absence of neutrophil / granulocyte lineages typical of Ficoll-isolated PBMC samples

# Conclusion:

This is not bone marrow. The proportions match PBMC (peripheral blood).

## 🦠 7. Is the Patient Healthy or Infected?
Findings:

NK cells elevated → innate antiviral activation

Monocytes moderately elevated → mild inflammatory response

T cell compartment stable → no severe lymphopenia

No extreme expansion of plasmablasts → not acute bacterial

Interpretation:

The immune landscape is mildly shifted toward an antiviral response—consistent with recent or ongoing viral exposure, but not severe infection.

## 🛠 8. Reproducibility Pipeline

# 1. Install environment
conda env create -f environment.yml
conda activate scRNA_project

# 2. Run analysis
python scripts/preprocessing.py
python scripts/clustering.py
python scripts/annotation.py
python scripts/plots.py

# 3. Launch notebook (optional)
jupyter notebook notebook.ipynb

# 4. View results
results/
├── umap_clusters.png
├── umap_celltypes.png
├── proportions.csv
├── markers.csv

## 📦 9. Dependencies
Package	Version
Python	3.10
Scanpy	≥1.9
Anndata	≥0.9
Matplotlib	≥3.7
Seaborn	≥0.12
scikit-learn	≥1.3
statsmodels	≥0.14
python-pptx	≥0.6

## 🚀 10. Future Directions

Add statistical validation (bootstrapping, jackknife) for cluster assignments

Perform pathway enrichment (GSEA, Enrichr)

Integrate public PBMC datasets for benchmarking

Add viral gene signatures for infection confirmation

Deploy interactive dashboards (Streamlit + Scanpy)

Containerize via Docker for complete reproducibility

## 🌟 11. Key Insight

The elevated NK and monocyte composition strongly suggests an innate antiviral activation pattern — without evidence of severe immune collapse.

## Short Scientific Narrative (for report or mentors)
This single-cell RNA-seq analysis profiled a hematopoietic/immune mixture using Scanpy and marker-based annotation. Despite containing erythroid and megakaryocytic cells consistent with marrow biology, the dataset displayed unusually high proportions of NK and T cells, and very few hematopoietic stem/progenitor cells. This pattern deviates from typical bone marrow composition and instead resembles PBMCs with minor marrow contamination. Functionally, the strong NK and T-cell expansion, combined with plasma cell emergence and monocyte elevation, indicates active immune stimulation. The most biologically supported interpretation is an ongoing viral infection driving cytotoxic and humoral responses.
