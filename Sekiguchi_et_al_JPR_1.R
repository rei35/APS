library(scDblFinder)
library(Seurat)

##Doublet/multiplet removal

samples <-c("BM150","BM152","BM156","BM157","BM158","BM165","BM168","BM169","GM136","GM143","GM144","GM147","GM148","GM169","GM183","GM184a","GM238","GM241","GM242","GM283","GM289")
for (element in samples) {
  set.seed(432)
  sample <- readRDS(paste0("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/",element,".RDS"))
  sample_sce <- as.SingleCellExperiment(sample)
  sample_sce <- scDblFinder(sample_sce)
  sample_Dbl <- as.Seurat(sample_sce, counts = "counts")
  Idents(sample_Dbl) <- "scDblFinder.class"
  sample_clean <- subset(sample_Dbl, idents = "singlet")
  saveRDS(sample_clean, file = paste0(element,"_ready.RDS"))
  rm(sample,sample_sce,sample_Dbl, sample_clean)
  gc()
}

##Load all samples and add sample metadata

BM150 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM150_ready.RDS")
BM150$sample <- rep("BM150",length(BM150$orig.ident))
BM150$tissue <- rep("Buccal",length(BM150$orig.ident))

BM152 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM152_ready.RDS")
BM152$sample <- rep("BM152",length(BM152$orig.ident))
BM152$tissue <- rep("Buccal",length(BM152$orig.ident))

BM156 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM156_ready.RDS")
BM156$sample <- rep("BM156",length(BM156$orig.ident))
BM156$tissue <- rep("Buccal",length(BM156$orig.ident))

BM157 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM157_ready.RDS")
BM157$sample <- rep("BM157",length(BM157$orig.ident))
BM157$tissue <- rep("Buccal",length(BM157$orig.ident))

BM158 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM158_ready.RDS")
BM158$sample <- rep("BM158",length(BM158$orig.ident))
BM158$tissue <- rep("Buccal",length(BM158$orig.ident))

BM165 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM165_ready.RDS")
BM165$sample <- rep("BM165",length(BM165$orig.ident))
BM165$tissue <- rep("Buccal",length(BM165$orig.ident))

BM168 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM168_ready.RDS")
BM168$sample <- rep("BM168",length(BM168$orig.ident))
BM168$tissue <- rep("Buccal",length(BM168$orig.ident))

BM169 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/BM169_ready.RDS")
BM169$sample <- rep("BM169",length(BM169$orig.ident))
BM169$tissue <- rep("Buccal",length(BM169$orig.ident))

GM136 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM136_ready.RDS")
GM136$sample <- rep("GM136",length(GM136$orig.ident))
GM136$tissue <- rep("Gingiva",length(GM136$orig.ident))

GM143 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM143_ready.RDS")
GM143$sample <- rep("GM143",length(GM143$orig.ident))
GM143$tissue <- rep("Gingiva",length(GM143$orig.ident))

GM144 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM144_ready.RDS")
GM144$tissue <- rep("Gingiva",length(GM144$orig.ident))
GM144$sample <- rep("GM144",length(GM144$orig.ident))

GM147 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM147_ready.RDS")
GM147$tissue <- rep("Gingiva",length(GM147$orig.ident))
GM147$sample <- rep("GM147",length(GM147$orig.ident))

GM148 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM148_ready.RDS")
GM148$tissue <- rep("Gingiva",length(GM148$orig.ident))
GM148$sample <- rep("GM148",length(GM148$orig.ident))

GM169 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM169_ready.RDS")
GM169$tissue <- rep("Gingiva",length(GM169$orig.ident))
GM169$sample <- rep("GM169",length(GM169$orig.ident))

GM183 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM183_ready.RDS")
GM183$tissue <- rep("Gingiva",length(GM183$orig.ident))
GM183$sample <- rep("GM183",length(GM183$orig.ident))

GM184a <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM184a_ready.RDS")
GM184a$tissue <- rep("Gingiva",length(GM184a$orig.ident))
GM184a$sample <- rep("GM184a",length(GM184a$orig.ident))

GM238 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM238_ready.RDS")
GM238$tissue <- rep("Gingiva",length(GM238$orig.ident))
GM238$sample <- rep("GM238",length(GM238$orig.ident))

GM241 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM241_ready.RDS")
GM241$tissue <- rep("Gingiva",length(GM241$orig.ident))
GM241$sample <- rep("GM241",length(GM241$orig.ident))

GM242 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM242_ready.RDS")
GM242$tissue <- rep("Gingiva",length(GM242$orig.ident))
GM242$sample <- rep("GM242",length(GM242$orig.ident))

GM283 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM283_ready.RDS")
GM283$tissue <- rep("Gingiva",length(GM283$orig.ident))
GM283$sample <- rep("GM283",length(GM283$orig.ident))

GM289 <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/GM289_ready.RDS")
GM289$sample <- rep("GM289",length(GM289$orig.ident))
GM289$tissue <- rep("Gingiva",length(GM289$orig.ident))

all <- merge(BM150, y= c(BM152,BM156,BM157,BM158,BM165,BM168,BM169,BM169,GM136,GM143,GM144,GM147,GM148,GM169,GM183,GM184a,GM238,GM241,GM242,GM283,GM289),
             add.cell.ids = c("BM150","BM152","BM156","BM157","BM158","BM165","BM168","BM169","BM169","GM136","GM143","GM144","GM147","GM148","GM169","GM183","GM184a","GM238","GM241","GM242","GM283","GM289"),
             project = "Oral_mucosa")
saveRDS(all,file = "all.RDS",compress = FALSE)

##Reading object with all samples

all <- readRDS("~/YamadaK_lab/Rei_2022/Mucosa/Dataset/all.RDS")
all

##Calculating mitochondrial DNA content and filtering cells out

all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^MT-")
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
all <- subset(all, subset = nFeature_RNA >= 200 & percent.mt <=12.5)

##Clearing some RAM

rm(BM150,BM152,BM156,BM157,BM158,BM165,BM168,BM169,BM169,GM136,GM143,GM144,GM147,GM148,GM169,GM183,GM184a,GM238,GM241,GM242,GM283,GM289)
gc()

##Normalization and clustering

all <- SCTransform(all, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
gc()
all <- RunPCA(all, verbose = TRUE)
ElbowPlot(all, ndims = 50)
all <- RunUMAP(all, dims = 1:50, verbose = FALSE)

all <- FindNeighbors(all, dims = 1:50)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{all <- FindClusters(all, resolution = resolution)
}

##Optimal clustering determination

library(clustree)

pdf('all_clustree_seurat.pdf', width = 30, height = 20)
clustree(all,node_size_range=c(10,20), node_text_size = 8)
dev.off()

##Selecting optimal resolution

Idents(object=all) <- "SCT_snn_res.0.1"

##Generating cluster proportions

mytable<-table(Idents(all),all$tissue)
mydf <- as.data.frame(prop.table(mytable,2))


pdf('cluster_proportions.pdf', width = 30, height = 20)
ggplot(mydf, aes(x=Var2, y = Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "Cluster") +
  xlab("tissue") + 
  ylab("Proportion")
dev.off()

##Visualizing clusters

pdf('Seurat_clusters.pdf', width = 30, height = 20)
DimPlot(all, reduction = "umap", label = TRUE, raster = FALSE, shuffle = FALSE)
DimPlot(all, reduction = "umap", group.by = "sample", raster = FALSE, shuffle = FALSE) + ggplot2::theme(legend.position = "bottom")
DimPlot(all, reduction = "umap", group.by = "tissue", raster = FALSE, shuffle = FALSE) + ggplot2::theme(legend.position = "bottom")
DimPlot(all, reduction = "umap", split.by = "sample", ncol = 3, raster = FALSE, shuffle = FALSE) + ggplot2::theme(legend.position = "bottom")
DimPlot(all, reduction = "umap", split.by = "tissue", ncol = 3, raster = FALSE, shuffle = FALSE) + ggplot2::theme(legend.position = "bottom")
DimPlot(all, reduction = "umap", group.by= "sample", split.by = "tissue", ncol = 3, label=TRUE, raster = FALSE, shuffle = FALSE) + ggplot2::theme(legend.position = "bottom")
dev.off()

table(all$orig.ident)
table(all$seurat_clusters)

##Visualizing markers

pdf('Seurat_cell_markers.pdf', width = 10, height = 20)
FeaturePlot(all, features = c("ACKR1","RAMP2","SELE","VWF","PECAM1","LUM","COL3A1","DCN","COL1A1","CFD","KRT14","KRT5","S100A2","CSTA","SPRR1B","CD69","CD52","CXCR4","PTPRC","HCST"), ncol = 3, slot = "counts", order = TRUE)
FeaturePlot(all, features = c("LUM","COL3A1","DCN","COL1A1"), ncol = 3, order = TRUE, raster = T) 
FeaturePlot(all,features = c("KRT14","KRT5","KRT19","SOX10","COL1A1","FN1"), ncol = 3, slot = "counts") 
FeaturePlot(all,features = c("ANGPT2","SELP","VWF","PECAM1")) #Endothelial cells
dev.off()

##Creating a small subset for heatmap representation

small.all<-subset(x = all, downsample = 100)


##Identifying cluster markers

all.markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(openxlsx)
write.xlsx(all.markers, file = 'all_markers.xlsx', rowNames=TRUE)

all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pdf('markers_heatmap.pdf', width = 10, height = 20)
DoHeatmap(small.all, features = top10$gene) + NoLegend()
dev.off()

##Running DE for fibroblast clusters from gingiva and mucosa 

DefaultAssay(all)<-"RNA"
fibroblast_DE_genes<-FindMarkers(all, ident.1 = 3, ident.2 = 4, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(fibroblast_DE_genes, file = 'fibroblast_DE_genes.xlsx', rowNames=TRUE)
DefaultAssay(all)<-"SCT"
FeaturePlot(all,features = c("CXCL13","TNXB"), ncol = 2)

##Creating artificial Ident to be able to compare cell by cluster-tissue

all$cluster.tissue <- paste(Idents(all), all$tissue, sep = "_")
all$cluster <- Idents(all)
Idents(all) <- "cluster.tissue"
Idents(all)
DefaultAssay(all)<-"RNA"
epithelial_DE_genes <- FindMarkers(all, ident.1 = "5_Gingiva", ident.2 = "5_Buccal", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(epithelial_DE_genes, file = 'epithelial_DE_genes.xlsx', rowNames=TRUE)

