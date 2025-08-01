# CEACAM_scdownstream

### GSE128169 _matrix.mtx.gz 转成 rds (nf-core/scdownstream 输入格式之一)
```R
Rscript mtx2rds.R
```
### 运行 nf-core/scdownstream 流程
##### 准备samplesheet.csv, executor.config
```bash
samplesheet.csv
sample,unfiltered
SC45NOR,data/GSM3666096_SC45NOR.rds
SC56NOR,data/GSM3666097_SC56NOR.rds
SC59NOR,data/GSM3666098_SC59NOR.rds
SC155NORLOW,data/GSM3666099_SC155NORLOW.rds
SC156NORUP,data/GSM3666100_SC156NORUP.rds
SC51SSCLOW,data/GSM3666101_SC51SSCLOW.rds
SC52SSCUP,data/GSM3666102_SC52SSCUP.rds
SC63SSCLOW,data/GSM3666103_SC63SSCLOW.rds
SC64SSCUP,data/GSM3666104_SC64SSCUP.rds
SC108SSCLOW,data/GSM3666105_SC108SSCLOW.rds
SC109SSCUP,data/GSM3666106_SC109SSCUP.rds
SC135SSCLOW,data/GSM3666107_SC135SSCLOW.rds

executor.config
params {
    celltypist_model = 'Human_PF_Lung'
}
```

##### 运行nf-core/scdownstream流程，我用的是singularity环境，根据自己服务器环境选择profile (e.g. conda, docker)
```bash
launch.sh
nextflow run nf-core/scdownstream -r dev --input ./samplesheet.csv --outdir ./results -profile singularity -c ./executor.config -resume
```

##### 得到 celltypist 'Human_PF_Lung' 注释的结果
```bash
results/finalized/merged.h5ad
```

##### 分类SSC/NOR，提取细胞label，提取CEACAM5、CEACAM6基因、作图
```bash
>>> print(adata.obs['disease_status'].value_counts())
disease_status
Disease    42496
Normal     23998

>>> print(adata.obs['celltypist:Human_PF_Lung'].value_counts())
celltypist:Human_PF_Lung
Macrophages                       31441
T Cells                            5978
Endothelial Cells                  5610
Ciliated                           5087
AT2                                3406
Monocytes                          2336
SCGB3A2+ SCGB1A1+                  1544
Myofibroblasts                     1532
NK Cells                           1292
SCGB3A2+                           1086
Smooth Muscle Cells                1067
Proliferating Macrophages           868
Mast Cells                          862
MUC5B+                              804
Basal                               638
Lymphatic Endothelial Cells         600
Differentiating Ciliated            427
cDCs                                386
B Cells                             298
AT1                                 271
Plasma Cells                        270
Fibroblasts                         177
Transitional AT2                    172
PLIN2+ Fibroblasts                  129
Proliferating Epithelial Cells       78
pDCs                                 38
Mesothelial Cells                    34
KRT5-/KRT17+                         31
Proliferating T Cells                26
MUC5AC+ High                          4
HAS1 High Fibroblasts                 2
```

##### 结果
```bash
CEACAM5_boxplot_with_stats.png
CEACAM6_boxplot_with_stats.png
```
