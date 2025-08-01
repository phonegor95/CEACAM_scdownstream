# ---- 1. libraries ------------------------------------------------------------
library(Seurat)      # Seurat v5 works too
library(Matrix)      # for sparse matrices
library(dplyr)       # handy for a few tweaks

# ---- 2. point to the directory that holds the files --------------------------
data_dir <- "/mnt/8TB-HDD/hongyongfeng/20250707_scdownstream/GSE128169"      # adjust!
setwd(data_dir)

# ---- 3. discover all matrix files & infer sample prefixes --------------------
mtx_files  <- list.files(pattern = "_matrix\\.mtx\\.gz$")
samples    <- sub("_matrix\\.mtx\\.gz$", "", mtx_files)

# sanity-check
print(samples)

sample_dirs <- file.path(data_dir, samples)

objs <- lapply(sample_dirs, function(d) {
  
  # read 10X‑style trio inside directory d
  counts <- Read10X(data.dir = d, gene.column = 2)   # 2 = gene symbol
  
  # build Seurat object
  obj <- CreateSeuratObject(counts, project = basename(d))
  
  # keep track of origin & ensure unique barcodes
  obj$sample <- basename(d)
  obj        <- RenameCells(obj, add.cell.id = basename(d))
  
  # write to disk  ➜  /path/to/rds/GSM3666096_SC45NOR.rds, etc.
  saveRDS(obj, file = file.path(d, paste0(basename(d), ".rds")))
  
  obj   # returned so objs[[i]] is still usable later
})
