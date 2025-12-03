# HackBio-Single-Cell-RNA-Seq-Stage-2
### **Single-Cell RNA-seq Analysis of Bone Marrow Dataset Using Scanpy**
This repository reproduces a complete scRNA-seq analysis pipeline using the *Scanpy* library on a modified bone marrow dataset (originally from CZI). The workflow includes preprocessing, normalization, clustering, marker-based annotation, and biological interpretation.

### **Workflow Overview**
1. Load `.h5ad` dataset
2. QC filtering: n_genes, percent mitochondrial, doublets
3. Normalization & log1p transformation
4. Highly variable gene selection
5. PCA → neighbors → UMAP
6. Leiden clustering
7. Marker gene scoring using **decoupler** with corrected Ensembl gene names
8. Cell-type annotation using PanglaoDB markers
9. Biological interpretation of tissue source & immune status
10. Visualization and reporting

### **Key Results**
- Identified major immune and hematopoietic populations.
- NK (36%) and T cells (32%) dominate the landscape.
- Presence of erythroid precursors & megakaryocytes suggests marrow origin.
- Relative absence of HSCs & lymphocyte-high profile suggests PBMC contamination.
- Immune activation signature consistent with acute viral infection.

### **How to Run**
```
!pip install scanpy decoupler-py omnipath pandas numpy matplotlib seaborn
google colab
```
Open the provided notebook and run cells sequentially.

### **Files in Repository**
- `analysis.ipynb` – Full Scanpy pipeline
- `cluster_markers_top10.csv` – Top marker genes per cluster
- `cluster_proportions.csv` – Proportion of annotated clusters
- `README.md` – Project overview & instructions

---

## Intro, Results, and Conclusion (Notebook Narrative)

### **Introduction**
Single-cell RNA-sequencing enables the characterization of heterogeneous cell populations at high resolution. In this project, we reproduced a standard Scanpy workflow using a modified bone marrow dataset. The goal was to identify major hematopoietic and immune cell types, evaluate the tissue origin of the sample, and infer the patient’s immune status from the cellular composition.

### **Results**
Dimensionality reduction and Leiden clustering revealed clear separation of NK cells, T cells, monocytes, B cells, plasma cells, erythroid precursors, megakaryocytes, and rare progenitors. Marker-based annotation using PanglaoDB (corrected for Ensembl IDs) confirmed these identities.

The composition showed striking NK (36%) and T-cell (32%) expansion, moderate monocyte elevation (13.8%), and a detectable plasma cell compartment (5.5%). Erythroid and megakaryocyte populations indicated marrow linkage, but extremely low HSC levels and lymphocyte dominance suggested a PBMC-like sample rather than full bone marrow.

### **Conclusion**
The dataset does not resemble typical bone marrow, which is usually myeloid-rich and progenitor-rich. Instead, it displays a lymphocyte-skewed profile indicative of an activated immune response. The NK/T-cell expansion, along with increased plasma cells, suggests ongoing **viral infection** rather than a healthy immune state. Overall, the analysis demonstrates how scRNA-seq data composition can reveal both tissue origin and real-time immunological events.

---

## Short Scientific Narrative (for report or mentors)
This single-cell RNA-seq analysis profiled a hematopoietic/immune mixture using Scanpy and marker-based annotation. Despite containing erythroid and megakaryocytic cells consistent with marrow biology, the dataset displayed unusually high proportions of NK and T cells, and very few hematopoietic stem/progenitor cells. This pattern deviates from typical bone marrow composition and instead resembles PBMCs with minor marrow contamination. Functionally, the strong NK and T-cell expansion, combined with plasma cell emergence and monocyte elevation, indicates active immune stimulation. The most biologically supported interpretation is an ongoing viral infection driving cytotoxic and humoral responses.
