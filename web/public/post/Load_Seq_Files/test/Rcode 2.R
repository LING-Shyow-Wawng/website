library(Seurat)
library(parallel)
library(future)
options(future.globals.maxSize =1024*1024^2)
library(stringr)
library(tidyverse)
#multiprocess(workers = 4) #don't set a number exceeding the maximum process

rm(integrate_objects)

save.image(file = "project_image.RData")

#Read each seq file
sham1 <- ReadMtx(
  mtx="../GSE174574_RAW/GSM5319987_sham1_matrix.mtx.gz",
  cells = "../GSE174574_RAW/GSM5319987_sham1_barcodes.tsv.gz",
  features = "../GSE174574_RAW/GSM5319987_sham1_genes.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
sham1_seurat <- CreateSeuratObject(counts = sham1)
sham1_seurat$condition <- "sham"
sham1_seurat[["percent_mt"]] <- PercentageFeatureSet(sham1_seurat, pattern = "^mt-")

sham2 <- ReadMtx(
  mtx="../GSE174574_RAW/GSM5319988_sham2_matrix.mtx.gz",
  cells = "../GSE174574_RAW/GSM5319988_sham2_barcodes.tsv.gz",
  features = "../GSE174574_RAW/GSM5319988_sham2_genes.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
sham2_seurat <- CreateSeuratObject(counts = sham2)
sham2_seurat$condition <- "sham"
sham2_seurat[["percent_mt"]] <- PercentageFeatureSet(sham2_seurat, pattern = "^mt-")

sham3 <- ReadMtx(
  mtx="../GSE174574_RAW/GSM5319989_sham3_matrix.mtx.gz",
  cells = "../GSE174574_RAW/GSM5319989_sham3_barcodes.tsv.gz",
  features = "../GSE174574_RAW/GSM5319989_sham3_genes.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
sham3_seurat <- CreateSeuratObject(counts = sham3)
sham3_seurat$condition <- "sham"
sham3_seurat[["percent_mt"]] <- PercentageFeatureSet(sham3_seurat, pattern = "^mt-")

stroke1 <- ReadMtx(
  mtx="../GSE174574_RAW/GSM5319990_MCAO1_matrix.mtx.gz",
  cells = "../GSE174574_RAW/GSM5319990_MCAO1_barcodes.tsv.gz",
  features = "../GSE174574_RAW/GSM5319990_MCAO1_genes.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
stroke1_seurat <- CreateSeuratObject(counts = stroke1)
stroke1_seurat$condition <- "stroke"
stroke1_seurat[["percent_mt"]] <- PercentageFeatureSet(stroke1_seurat, pattern = "^mt-")

stroke2 <- ReadMtx(
  mtx="../GSE174574_RAW/GSM5319991_MCAO2_matrix.mtx.gz",
  cells = "../GSE174574_RAW/GSM5319991_MCAO2_barcodes.tsv.gz",
  features = "../GSE174574_RAW/GSM5319990_MCAO1_genes.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
stroke2_seurat <- CreateSeuratObject(counts = stroke2)
stroke2_seurat$condition <- "stroke"
stroke2_seurat[["percent_mt"]] <- PercentageFeatureSet(stroke2_seurat, pattern = "^mt-")

stroke3 <- ReadMtx(
  mtx="../GSE174574_RAW/GSM5319992_MCAO3_matrix.mtx.gz",
  cells = "../GSE174574_RAW/GSM5319992_MCAO3_barcodes.tsv.gz",
  features = "../GSE174574_RAW/GSM5319992_MCAO3_genes.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
stroke3_seurat <- CreateSeuratObject(counts = stroke3)
stroke3_seurat$condition <- "stroke"
stroke3_seurat[["percent_mt"]] <- PercentageFeatureSet(stroke3_seurat, pattern = "^mt-")

#Create a list to store all seurat files
samples_objects <- c(sham1_seurat, sham2_seurat, sham3_seurat, stroke1_seurat, stroke2_seurat, stroke3_seurat)
sample_names <- c("sham1","sham2","sham3","stroke1","stroke2","stroke3")
names(samples_objects) <- sample_names

# Analysis parameters
mc <- 4                    # number of given cores
anchor_dims <- 30            # number of anchor dimensions used for biological replicates integration
pca_dims <- 30               # number of PCA dimensions to compute and use in tSNE, UMAP 
umap_n_neighbors <- 30       # UMAP parameter
clustering_resolution <- 0.3 # Resolution parameter for Seurat clustering
n_features <- 2000

# Filter cells (based on percent_mito and nFeature_RNA), normalize (log-normalize in columns) 
# and find n_features most variable features/genes (using Variance-stabilizing transformation)
# to the list of top variable genes add markers of microglia, macrophages and cell cycle
samples_objects <- mclapply(seq_along(samples_objects), function(i) {
  samples_objects[[i]] <- subset(x = samples_objects[[i]],
                                 subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt < 5)
  samples_objects[[i]] <- NormalizeData(object = samples_objects[[i]], verbose = FALSE, scale.factor = 1e6)
  samples_objects[[i]] <- FindVariableFeatures(object = samples_objects[[i]],
                                               selection.method = "vst", 
                                               nfeatures = n_features, verbose = FALSE)
  # Add selected (previously reported) genes to var.features
  samples_objects[[i]]@assays$RNA@var.features <- unique(
    c(samples_objects[[i]]@assays$RNA@var.features, 
      s_genes = str_to_title(cc.genes.updated.2019$s.genes),
      g2m_genes= str_to_title(cc.genes.updated.2019$g2m.genes)))
  samples_objects[[i]]
}, mc.cores = "3")
names(samples_objects) <- sample_names

rm(integrate_objects)
# Integrate replicates within conditions, scale, regress out unwanted sources of variation, 
# calculate PCA, t-SNE and find cell clusters
integrate_objects <- mclapply(c(1), function(i) { #sham1, stroke1
  print("Find Integration Anchors: ")
  samples_anchors <- FindIntegrationAnchors(object.list = list(samples_objects[[i]], 
                                                               samples_objects[[i + 1]],
                                                               samples_objects[[i + 2]],
                                                               samples_objects[[i + 3]],
                                                               samples_objects[[i + 4]],
                                                               samples_objects[[i + 5]]), 
                                            dims = 1:anchor_dims,
                                            anchor.features = n_features)
  print("Integrate Data: ")
  samples_integrated <- IntegrateData(anchorset = samples_anchors, dims = 1:anchor_dims)
  DefaultAssay(object=samples_integrated) <- "integrated"
  
  print("Cell Cycle Scoring:")
  samples_integrated <- CellCycleScoring(object = samples_integrated, 
                                         s.features = str_to_title(cc.genes.updated.2019$s.genes), 
                                         g2m.features = str_to_title(cc.genes.updated.2019$g2m.genes), 
                                         set.ident = TRUE)  
  
  print("Scale data, regress out: nCount_RNA, percent_mt, CC_difference:")
  
  # Calculate difference between G2M and S phase scores to separate non-cycling and cycling cells
  # Approch described in Seurat's Cell-Cycle Scoring and Regression vignette 
  samples_integrated$CC_Difference <- samples_integrated$S.Score - samples_integrated$G2M.Score
  samples_integrated <- ScaleData(object = samples_integrated, verbose = FALSE, 
                                  vars.to.regress = c("nCount_RNA", 
                                                      "percent_mt", 
                                                      "CC_Difference"))  
  
  print("Run PCA: ")
  samples_integrated <- RunPCA(object = samples_integrated, npcs = pca_dims, verbose = FALSE)
  
  print("Run t-SNE:")
  samples_integrated <- RunTSNE(object = samples_integrated, reduction = "pca", dims=1:pca_dims)
  
  print("Run UMAP:")
  samples_integrated <- RunUMAP(object = samples_integrated, reduction = "pca", dims=1:pca_dims)
  
  print("Find Neighbors: ")
  samples_integrated <- FindNeighbors(object = samples_integrated, dims = 1:pca_dims)
  
  print("Find Clusters: ")
  samples_integrated <- FindClusters(object = samples_integrated, resolution = clustering_resolution) 
  samples_integrated
}, mc.cores = "3")

rm(cluster_markers)
# Find markers for obtained clusters
cluster_markers <- FindAllMarkers(integrated_objects, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 1)
markers16 <- FindMarkers(object = integrated_objects, 
                         ident.1 = 16,
                         only.pos = TRUE, 
                         min.pct = 0.25, 
                         logfc.threshold = 1)
markers17 <- FindMarkers(object = integrated_objects, 
                         ident.1 = 17,
                         only.pos = TRUE, 
                         min.pct = 0.25, 
                         logfc.threshold = 1)

DefaultAssay(integrated_objects) <- "RNA"

Idents(integrated_objects) <- "clustertype"
DimPlot(integrated_objects, reduction = "umap", group.by = "condition")
DimPlot(integrated_objects, reduction = "umap", label = T, repel = T)
DimPlot(integrated_objects, reduction = "umap", split.by = "condition", label = T, repel = T)

FeaturePlot(integrated_objects, 
            features = c("Tmem119", #microglia
                         "Hexb", #microglia
                         "Ccr2", #monocyte
                         "Mrc1", #macrophages
                         "Gfap", #astrocyte
                         "Aldh1l1", #astrocyte
                         "Olig2", #oligo
                         "Plp1", #oligo
                         "Ttr", #ependymocytes
                         "Itm2a",#endothelial cell
                         "Acta2",  #vascular smooth muscle cells (SMC)
                         "Gzmb", "Nkg7", #NKT
                         "Ttr", "Clu", #NK
                         "S100a8", "S100a9", #DCs,
                         "Col1a1", #Fibroblast
                         "Cd19", #Bcell
                         "Cd8a", "Cd8b1" ),#Tcell 
            cols = c("grey", "red"),
            label = F)


#Plot C3ar1 gene expression
FeaturePlot(integrated_objects, 
            features = "C3ar1",
            cols = c("grey", "red"),
            label = F)
FeaturePlot(integrated_objects, 
            features = "C3ar1",
            split.by = "condition",
            cols = c("grey", "red"),
            label = F)

#Plot C3 gene expression
FeaturePlot(integrated_objects, 
            features = "C3",
            cols = c("grey", "red"),
            label = F)
FeaturePlot(integrated_objects, 
            features = "C3",
            split.by = "condition",
            cols = c("grey", "red"),
            label = F)
VlnPlot(integrated_objects, 
        features = "C3", 
        split.by = "condition",
        pt.size = 0.0001, 
        combine = FALSE)

#Plot Cr2 gene expression
FeaturePlot(integrated_objects, 
            features = "Cr2",
            cols = c("grey", "red"),
            label = F)
FeaturePlot(integrated_objects, 
            features = "Cr2",
            split.by = "condition",
            cols = c("grey", "red"),
            label = F)
VlnPlot(integrated_objects, 
        features = "Cr2", 
        split.by = "condition",
        pt.size = 0.0001, 
        combine = FALSE)

#Plot Cd18 gene expression
FeaturePlot(integrated_objects, 
            features = "Itgb2",
            cols = c("grey", "red"),
            label = F)
FeaturePlot(integrated_objects, 
            features = "Itgb2",
            split.by = "condition",
            cols = c("grey", "red"),
            label = F)
VlnPlot(integrated_objects, 
        features = "Itgb2", 
        split.by = "condition",
        pt.size = 0.0001, 
        combine = FALSE)

#Plot Cd11b gene expression
FeaturePlot(integrated_objects, 
            features = "Itgam",
            cols = c("grey", "red"),
            label = F)
FeaturePlot(integrated_objects, 
            features = "Itgam",
            split.by = "condition",
            cols = c("grey", "red"),
            label = F)
VlnPlot(integrated_objects, 
        features = "Itgam", 
        split.by = "condition",
        pt.size = 0.0001, 
        combine = FALSE)

#plot cd11b/cd18 blend
FeaturePlot(integrated_objects, 
            features = c("Itgam","Itgb2"),
            blend = T,
            blend.threshold = 0.5,pt.size = 1,
            split.by = "condition",
            cols = c("grey", "light green", "dark red"),
            label = F)


Idents(integrate_objects) <- "clusterid"
# Rename classes.
integrated_objects <- RenameIdents(object = integrated_objects, 
                                   `0` = "Endothelia", 
                                   `1` = "Microglia", 
                                   `2` = "Endothelia",
                                   `3` = "Microglia",
                                   `4` = "Endothelia",
                                   `5` = "Monocytes",
                                   `6` = "Astrocytes",
                                   `7` = "CNS-associated Macrophages",
                                   `8` = "Natural Killer Cells",
                                   `9` = "Smooth Muscle Cells",
                                   `10` = "Pericytes",
                                   `11` = "Oligodendrocytes",
                                   `12` = "Dendritic Cells",
                                   `13` = "Microglia",
                                   `14` = "Lymphocytes",
                                   `15` = "Fibroblasts",
                                   `16` = "Capillary Endothelia",
                                   `17` = "Ependymal Cells",
                                   `18` = "Red Blood Cells")
Idents(integrated_objects)
FeaturePlot(integrated_objects, 
            features = "C3ar1",
            split.by = "condition",
            cols = c("grey", "red"),
            label = F)

VlnPlot(integrated_objects, 
        features = "C3ar1", 
        split.by = "condition",
        pt.size = 0.0001, 
        combine = FALSE)

#calculate DEGs of same cluster between sham and stroke conditions
#1.Store another metadata of cluster_id+ctrl/tumor
integrated_objects$clustertype_condition <- paste(Idents(integrated_objects), integrated_objects$condition, sep = "_")
#2.Store clusterid
integrated_objects$clustertype <- Idents(integrated_objects)
#3.Change cluster idents to clusterid_condition
Idents(integrated_objects) <- "clustertype_condition"
Idents(integrated_objects)
#4.Calcuate DEGs
DefaultAssay(integrated_objects) <- "RNA"
unique(integrated_objects$clustertype)
DEGs_by_cluster <- list() #store DEGs within clusters in a list 
for (i in 1:15) {
  group1 <- as.character(paste(unique(integrated_objects$clustertype)[i], "stroke", sep = "_"))
  group2 <- as.character(paste(unique(integrated_objects$clustertype)[i], "sham", sep = "_"))
  DEGs_by_cluster[[i]] <- FindMarkers(integrated_objects, 
                                      ident.1 = group1, 
                                      ident.2 = group2,
                                      only.pos = T,
                                      min.pct = 0.25, 
                                      verbose = FALSE)
}

monocyte_degs <- FindMarkers(integrated_objects, 
                             ident.1 = as.character(paste(unique(integrated_objects$clustertype)[2], "stroke", sep = "_")),
                             ident.2 = as.character(paste(unique(integrated_objects$clustertype)[2], "sham", sep = "_")),
                             min.pct = 0.25, 
                             verbose = FALSE)

isC3ar1 <- list() #inquire DEGs within clusters in a list 
for (i in c(1:15)) {
  isC3ar1[[i]] <- grep(rownames(DEGs_by_cluster[[i]]), pattern = "C3ar1")
} #only microglia & CNS-Macrophages significantly upregulated 
DEGs_by_cluster[[4]]["C3ar1",] #microglia
DEGs_by_cluster[[6]]["C3ar1",] #CNS-macrophages
DEGs_by_cluster[[15]]["C3ar1",] #CNS-macrophages
DEGs_by_cluster[[2]]["C3ar1",] #microglia

#cell velocity
library(devtools)
options(buildtools.check = function(action) TRUE )
BiocManager::install("pcaMethods") #velocity.R dependents
install_github("velocyto-team/velocyto.R")
remotes::install_github('satijalab/seurat-wrappers')

#Pseudotime 
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
Idents(integrated_objects) <- "clustertype"
#extract microglia and MoM cluster
micro_seurat <- subset(integrated_objects,
                       idents = "Microglia")
micro_seurat <- FindVariableFeatures(object = micro_seurat,
                                             selection.method = "vst", 
                                             nfeatures = n_features, verbose = FALSE)
micro_seurat <- ScaleData(micro_seurat)
micro_seurat <- RunPCA(object = micro_seurat, npcs = pca_dims, verbose = FALSE)
micro_seurat <- RunUMAP(object = micro_seurat, reduction = "pca", dims=1:3)
micro_seurat <- FindNeighbors(object = micro_seurat, dims = 1:3)
micro_seurat <- FindClusters(object = micro_seurat, resolution = 0.1) 
DimPlot(micro_seurat)
DimPlot(micro_seurat, group.by = "condition")
FeaturePlot(micro_seurat, "C3ar1", cols = c("grey", "red"))

microglia.cds <- as.cell_data_set(micro_seurat)
microglia.cds <- cluster_cells(cds = microglia.cds, reduction_method = "UMAP")
microglia.cds <- learn_graph(microglia.cds, use_partition = TRUE)
microglia.cds <- order_cells(microglia.cds, reduction_method = "UMAP")
plot_cells(
  cds = microglia.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  graph_label_size = 5
)
mom_seurat <- subset(integrated_objects,
                     idents = c("Monocytes","CNS-associated Macrophages"))
