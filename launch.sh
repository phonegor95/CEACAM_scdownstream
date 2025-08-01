nextflow run nf-core/scdownstream -r dev --input ./samplesheet.csv --outdir ./results -profile singularity -c ./executor.config -resume

#nextflow run nf-core/scdownstream -r dev --base_adata results/finalized/merged_cleaned_filtered.h5ad --base_embeddings scvi --base_label_col celltypist_clean --cluster_per_label true --outdir ./results -profile singularity -c ./executor.config -resume
