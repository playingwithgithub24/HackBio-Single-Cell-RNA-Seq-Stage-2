# A simple ISG signature check (common ISGs)
ISG_genes = ['IFIT1','IFIT2','IFIT3','MX1','ISG15','OAS1','OAS2','RSAD2','IFI6','IFI27','IFI44L']

# Create a gene name to Ensembl ID mapping dictionary (if not already available, assume ensembl_var is global)
gene_name_to_ensembl = pd.Series(ensembl_var.ensembl_gene_id.values, index=ensembl_var.gene_name).to_dict()

# Convert ISG gene symbols to Ensembl IDs
isg_ensembl_ids = []
for gene_symbol in ISG_genes:
    ensembl_id = gene_name_to_ensembl.get(gene_symbol)
    if ensembl_id: # Only add if a mapping was found
        isg_ensembl_ids.append(ensembl_id)

# Filter for ISG genes present in adata.var_names (now using Ensembl IDs)
isg_present = [g for g in isg_ensembl_ids if g in adata.var_names]
print("ISG genes present (Ensembl IDs):", isg_present)

# Compute module score for ISG signature per cell and compare across annotated cell types
# Only run if there are valid genes to score
if len(isg_present) > 0:
    sc.tl.score_genes(adata, gene_list=isg_present, score_name='ISG_score')
    # Save average ISG score per cell type
    isg_by_celltype = adata.obs.groupby('celltype_suggested')['ISG_score'].mean()
    write_csv(isg_by_celltype.rename("ISG_score").to_frame(), "ISG_by_celltype.csv")
    display(isg_by_celltype)
else:
    print("No ISG genes found in adata.var_names after Ensembl ID mapping. Skipping scoring.")
