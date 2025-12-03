# Single-Cell Analysis Report
Dataset: bone_marrow.h5ad

## 1. Identified cell types (annotated clusters)
- Cluster 0: T cells
- Cluster 1: NK cells
- Cluster 10: Erythroid
- Cluster 11: Plasma cells
- Cluster 12: T cells
- Cluster 13: Plasma cells
- Cluster 14: NK cells
- Cluster 15: Monocytes
- Cluster 16: Megakaryocytes/Platelets
- Cluster 17: NK cells
- Cluster 18: HSC/Progenitor
- Cluster 2: NK cells
- Cluster 3: B cells
- Cluster 4: Monocytes
- Cluster 5: T cells
- Cluster 6: Monocytes
- Cluster 7: NK cells
- Cluster 8: T cells
- Cluster 9: NK cells

## 2. Biological role of each cell type (short)
- T cells: Adaptive immune cells mediating cellular immunity (CD8 cytotoxic, CD4 helper).)
- NK cells: Innate lymphoid cells that kill virally infected or transformed cells rapidly without prior sensitization.)
- B cells: Adaptive immune cells producing antibodies (when differentiated to plasma cells) and presenting antigen.)
- Monocytes: Circulating precursors that can migrate to tissues and differentiate into macrophages or DCs; important for phagocytosis and inflammation.)
- Erythroid: Red blood cell precursors producing hemoglobin; present in marrow as erythroblasts.)
- Plasma cells: Terminally differentiated B cells that secrete large amounts of antibodies.)
- Megakaryocytes/Platelets: Megakaryocytes are large marrow cells that produce platelets — fragments that support clotting.)
- HSC/Progenitor: Hematopoietic stem/progenitor cells: give rise to all blood lineages and are typically found in bone marrow.)

## 3. Is the tissue source really bone marrow? Reasoning:
Provide an evidence-driven argument using expected vs missing lineage populations, presence of progenitors, and typical frequency distributions. Use the cluster proportions table in results/cluster_proportions.csv to support claims.

## 4. Health assessment (based on abundances)
Use neutrophil, monocyte, NK, and lymphocyte proportions to argue whether the donor likely has infection, inflammation, or is healthy. Provide quantitative thresholds (e.g., neutrophil overrepresentation > expected marrow proportion suggests emergency granulopoiesis).