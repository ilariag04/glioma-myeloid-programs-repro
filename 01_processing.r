library(dplyr)
library(Seurat)
library(sctransform)

# 1. Configurazione memoria
options(future.globals.maxSize = 64 * 1024^3)
options(bitmapType = 'cairo')
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

# 2. Caricamento dati
data_dir <- "/home/jovyan/work/shared/SCP2389/expression/66f44f68d5b3cc83dd766a0f"
counts <- Read10X(data.dir = data_dir, gene.column = 1)
Tumors.combined <- CreateSeuratObject(counts = counts, project = "SeqWell", min.cells = 3, min.features = 200)

# Libera subito la RAM della matrice grezza
rm(counts)
gc()

Tumors.combined[["percent.mt"]] <- PercentageFeatureSet(Tumors.combined, pattern = "^MT-")

# 3. Filtraggio QC
Tumors.combined <- subset(Tumors.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.mt < 20)

pdf("SeqWell_WT_Mutant_Tumors_QC_AF.pdf", height = 6, width = 20)
VlnPlot(Tumors.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# 4. Normalizzazione SCTransform e PCA
Tumors.combined <- SCTransform(Tumors.combined, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = TRUE)
Tumors.combined <- RunPCA(Tumors.combined)

pdf("SeqWell_WT_Mutant_Tumors_ElbowPlot.pdf", height = 6, width = 6)
ElbowPlot(Tumors.combined, ndims = 50)
dev.off()

write.table(Tumors.combined@meta.data, file = "SeqWell_WT_Mutant_Tumors_Integrated_MetaData.txt", sep = "\t", col.names = NA, quote = FALSE)

# 5. UMAP e Louvain Clustering
Tumors.combined <- RunUMAP(Tumors.combined, reduction = "pca", dims = 1:24)
Tumors.combined <- FindNeighbors(Tumors.combined, dims = 1:24)
Tumors.combined <- FindClusters(Tumors.combined, resolution = 0.3)

pdf("SeqWell_WT_Mutant_Tumors_Clusters_With_Labels.pdf", width = 10, height = 8)
DimPlot(Tumors.combined, reduction = "umap", label = TRUE)
dev.off()

# 6. Salvataggio Anti-OOM
Tumors.combined@assays$SCT@scale.data <- matrix(nrow = 0, ncol = 0)
gc()
saveRDS(Tumors.combined, file = "SeqWell_Brain_Tumors.rds", compress = FALSE)