library(scDblFinder)
library(Seurat)

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

Idents(all) <- "tissue"
table(Idents(all))

oral<-subset(all, idents = c("Buccal","Gingiva","Perio"))

Idents(oral) <- "celltype"

table(Idents(oral))

count <-0
levels(Idents(oral))
table(oral$celltype) ###sanity check
for (ident in levels(Idents(oral))) {
  if (length( WhichCells(oral, ident = ident))==1) {
    tempmatrix<-as.matrix(GetAssayData(oral, slot = "counts")[, WhichCells(oral, ident = ident)])
  } else {
    tempmatrix<-as.matrix(Matrix::rowSums(GetAssayData(oral, slot = "counts")[, WhichCells(oral, ident = ident)]))
  }
  count <- count + 1
  if (match(ident,levels(Idents(oral)))==1) {
    finalmatrix<-tempmatrix
  } else {
    finalmatrix<-cbind(finalmatrix,tempmatrix)
  }
}
colnames(finalmatrix)<-as.character(levels(Idents(oral)))
write.table(finalmatrix, file='~/Mucosa/Dataset1/oral_Raw_Gene_Counts_per_Celltype.tsv', quote=TRUE, sep='\t', col.names = TRUE)
val1=as.numeric(5) ## MINCOUNTS
val2=as.numeric(1) ## MINSAMPLES
filter <- apply(finalmatrix, 1, function(x) length(x[x>val1])>=val2)
fmtxfiltered=finalmatrix[filter,]
write.table(as.data.frame(fmtxfiltered),file='~/Mucosa/Dataset1/oral_Raw_Gene_Counts_per_Celltype_filtered.tsv', quote=TRUE, sep='\t', col.names = TRUE)


setwd("~/YamadaK_lab/Rei_2022/Mucosa/Dataset1")
set.seed(42)
library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
plan("multisession", workers = 4)
options(future.globals.maxSize = 6000 * 1024^2)



DefaultAssay(FibroblastsP)
FeaturePlot(FibroblastsP,features = c("CXCL9","CXCL10"), split.by = "tissue", ncol= 7, order = TRUE)

Idents(FibroblastsP)<-"tissue"
Idents(FibroblastsP)


FibroblastO <- readRDS("~/Mucosa/Dataset1/FibroblastO.RDS")

FibroblastO <- subset(FibroblastsP, idents = c("Buccal","Gingiva","Perio"))


DimPlot(FibroblastO, group.by = "tissue")
DimPlot(oral, group.by = "cell.type")

FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10"), split.by = "tissue", ncol= 7, order = TRUE)

FibroblastO <- RunPCA(FibroblastO, verbose = TRUE)
ElbowPlot(FibroblastO, ndims = 50)
FibroblastO <- RunUMAP(FibroblastO, dims = 1:50, verbose = FALSE)

FibroblastO <- FindNeighbors(FibroblastO, dims = 1:50)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{FibroblastO <- FindClusters(FibroblastO, resolution = resolution)
}

library(clustree)

pdf('FibroblastO_clustree_seurat.pdf', width = 30, height = 20)
clustree(FibroblastO,node_size_range=c(10,20), node_text_size = 8)
dev.off()

Idents(object=FibroblastO) <- "SCT_snn_res.0.1"

DimPlot(FibroblastO, group.by = "SCT_snn_res.0.3")

FeaturePlot(all,features = c("OCT4","SSEA4","STRO1","SOX2","TP64","TP75"), split.by = "tissue", ncol= 7, order = TRUE)


pdf('PERICYTE.pdf', width = 19, height = 5)
FeaturePlot(all,features = c("PDGFRB","MCAM","CSPG4","KCNJ8","ABCC9","MYH11","HIGD1B","KCNA5","ACTA2","CSPG4","PLN","RERGL","CD248"), ncol= 6, order = TRUE)
dev.off()


FeaturePlot(allP,features = c("PDGFRB","MCAM","CSPG4","KCNJ8","ABCC9","MYH11","HIGD1B","KCNA5"), ncol= 5, order = TRUE)
#OTHER PERICYTE MARKERS: "KCNK3",

FeaturePlot(FibroblastO,features = c("KRT1","KRT2","KRT13","SLPI","KRT17","ANXA1","MYL6","ERP29"), split.by = "tissue", ncol= 7, order = TRUE)
FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), split.by = "tissue", ncol= 7, order = TRUE)
FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), ncol= 2, order = TRUE)
FeaturePlot(FibroblastO,features = c("WNT5A"), ncol= 3, order = TRUE,split.by = "tissue")
FeaturePlot(FibroblastO,features = c("PITX1"), ncol= 3, order = TRUE,split.by = "tissue")


Idents(object=FibroblastO) <- "SCT_snn_res.0.3"

saveRDS(FibroblastO, file = "FibroblastO.RDS")

FibroblastO.markers_0.3_res <- FindAllMarkers(FibroblastO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(openxlsx)
write.xlsx(FibroblastO.markers_0.3_res, file = 'FibroblastO.markers_0.3_res.xlsx', rowNames=TRUE)


FibroblastO.markers_0.3_res %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pdf('FibroblastO_markers_heatmap_top10.pdf', width = 8, height = 10)
DoHeatmap(FibroblastO, features = top10$gene) + NoLegend()
dev.off()

VlnPlot(FibroblastO, features = c("CXCL9","CXCL10","PITX1","WNT5A"), group.by = "SCT_snn_res.0.3", split.by = "tissue")
VlnPlot(FibroblastO, features = c("CXCL10"), group.by = "SCT_snn_res.0.3", split.by = "tissue")

FibroblastO$CXCL <- "no"

cxcl9pos <- WhichCells(FibroblastO, expression = CXCL9 > 0.5)
cxcl10pos <- WhichCells(FibroblastO, expression = CXCL10 > 0.5)


Idents(FibroblastO) <- "CXCL"

FibroblastO <- SetIdent(FibroblastO, cells =c(cxcl9pos,cxcl10pos), value = 'yes')
Idents(FibroblastO, cells =cxcl10pos) <- "yes"
FibroblastO$CXCL <- Idents(FibroblastO)

FibroblastO.CXCL9.10.markers <- FindMarkers(FibroblastO, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(FibroblastO.CXCL9.10.markers, file = 'FibroblastO.CXCL9.10.DEG.xlsx', rowNames=TRUE)

##Subset of buccal mucosa

Idents(FibroblastO)<-"tissue"
Buccal_fibro<-subset(FibroblastO, idents = "Buccal")
Gingiva_fibro<-subset(FibroblastO, idents = "Gingiva")
Perio_fibro<-subset(FibroblastO, idents = "Perio")

Buccal_fibro$wnt5a <- "no"
Buccal_fibro$wnt5a

wnt5a1pos <- WhichCells(Buccal_fibro, expression = WNT5A > 0.5) ##1700 cells

Buccal_fibro ##total number of cells 7530
Idents(Buccal_fibro) <- "wnt5a"
Idents(Buccal_fibro)
Buccal_fibro <- SetIdent(Buccal_fibro, cells =wnt5a1pos, value = 'yes')
table(Idents(Buccal_fibro))
Buccal_fibro$wnt5a <- Idents(Buccal_fibro)
table(Idents(Buccal_fibro))

Buccal_fibro.wnt5a.markers <- FindMarkers(Buccal_fibro, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(Buccal_fibro.wnt5a.markers, file = 'Buccal_fibro.wnt5a.DEG.xlsx', rowNames=TRUE)
saveRDS(Buccal_fibro, file = "Buccal_fibro.RDS")

##Subset of Gingiva mucosa

Idents(FibroblastO)<-"tissue"
Idents(FibroblastO)
Gingiva_fibro<-subset(FibroblastO, idents = "Gingiva")
Gingiva_fibro$wnt5a <- "no"
Gingiva_fibro$wnt5a
cxcl9buc <- WhichCells(Buccal_fibro, expression = CXCL9 > 0.5) ##5838 cells
wnt5a1pos <- WhichCells(Gingiva_fibro, expression = WNT5A > 0.5) ##5838 cells
Gingiva_fibro ##total number of cells 7530
Idents(Buccal_fibro) <- "wnt5a"
Idents(Buccal_fibro)
Buccal_fibro <- SetIdent(Buccal_fibro, cells =wnt5a1pos, value = 'yes')
table(Idents(Buccal_fibro))

Buccal_fibro$wnt5a <- Idents(Buccal_fibro)
table(Idents(Buccal_fibro))

Buccal_fibro.wnt5a.markers <- FindMarkers(Buccal_fibro, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(Buccal_fibro.wnt5a.markers, file = 'Buccal_fibro.wnt5a.DEG.xlsx', rowNames=TRUE)
saveRDS(Buccal_fibro, file = "Buccal_fibro.RDS")

##Subset of Perio mucosa

Idents(FibroblastO)<-"tissue"
Idents(FibroblastO)
Perio_fibro<-subset(FibroblastO, idents = "Perio")
Perio_fibro$wnt5a <- "no"
Perio_fibro$wnt5a
wnt5a1pos <- WhichCells(Perio_fibro, expression = WNT5A > 0.5) ##3083 cells
Perio_fibro ##total number of cells 3978

Idents(Perio_fibro) <- "wnt5a"
Idents(Perio_fibro)
Perio_fibro <- SetIdent(Perio_fibro, cells =wnt5a1pos, value = 'yes')
table(Idents(Perio_fibro))
Perio_fibro$wnt5a <- Idents(Perio_fibro)
table(Idents(Perio_fibro))

Perio_fibro.wnt5a.markers <- FindMarkers(Perio_fibro, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(Perio_fibro.wnt5a.markers, file = 'Perio_fibro.wnt5a.DEG.xlsx', rowNames=TRUE)
saveRDS(Perio_fibro, file = "Perio_fibro.RDS")

##Perio_fibro has 3978 cells
DefaultAssay(Perio_fibro)
VlnPlot(Perio_fibro, features = "WNT5A", slot = "data")
VlnPlot(Perio_fibro, features = "PITX1", slot = "data")
VlnPlot(Perio_fibro, features = "CXCL9", slot = "data")
VlnPlot(Perio_fibro, features = "CXCL10", slot = "data")
VlnPlot(Perio_fibro, features = "WNT5A", slot = "data")

wnt5apos <- WhichCells(Perio_fibro, expression = WNT5A > 0.1, slot = "data") # 3083 cells
pitx1pos <- WhichCells(Perio_fibro, expression = PITX1 > 0.1) # 95 cells
cxcl9pos <- WhichCells(Perio_fibro, expression = CXCL9 > 0.1) # 263 cell
cxcl10pos <- WhichCells(Perio_fibro, expression = CXCL10 > 0.1) # 132 cells
cxcl9.10pos <- WhichCells(Perio_fibro, expression = CXCL9 > 0.1 & CXCL10 > 0.1) #54 cells

Gingiva_fibro
##Gingiva_fibro has 10358 cells
wnt5apos <- WhichCells(Gingiva_fibro, expression = WNT5A > 0.1) # 5838 cells
pitx1pos <- WhichCells(Gingiva_fibro, expression = PITX1 > 0.1) # 220 cells
cxcl9pos <- WhichCells(Gingiva_fibro, expression = CXCL9 > 0.1) # 150 cell
cxcl10pos <- WhichCells(Gingiva_fibro, expression = CXCL10 > 0.1) # 152 cells
cxcl9.10pos <- WhichCells(Gingiva_fibro, expression = CXCL9 > 0.1 & CXCL10 > 0.1) #13 cells

Perio_fibro
##Perio_fibro has 3978 cells
wnt5apos <- WhichCells(Perio_fibro, expression = WNT5A > 0.1) # 3083 cells
pitx1pos <- WhichCells(Perio_fibro, expression = PITX1 > 0.1) # 95 cells
cxcl9pos <- WhichCells(Perio_fibro, expression = CXCL9 > 0.1) # 263 cell
cxcl10pos <- WhichCells(Perio_fibro, expression = CXCL10 > 0.1) # 132 cells
cxcl9.10pos <- WhichCells(Perio_fibro, expression = CXCL9 > 0.1 & CXCL10 > 0.1) #54 cells

setwd("~/YamadaK_lab/Rei_2022/Mucosa/Dataset1")
set.seed(42)
library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
plan("multisession", workers = 4)
options(future.globals.maxSize = 6000 * 1024^2)

DefaultAssay(FibroblastsP)
FeaturePlot(FibroblastsP,features = c("CXCL9","CXCL10"), split.by = "tissue", ncol= 7, order = TRUE)

Idents(FibroblastsP)<-"tissue"
Idents(FibroblastsP)

FibroblastO <- readRDS("~/Mucosa/Dataset1/FibroblastO.RDS")
FibroblastO <- subset(FibroblastsP, idents = c("Buccal","Gingiva","Perio"))

DimPlot(FibroblastO, group.by = "tissue")

#Fibrosis markers
FeaturePlot(FibroblastO,features = c("CCN2","POSTN","LOXL2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral,features = c("POSTN","LOXL2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(FibroblastO,features = c("POSTN","LOXL2"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(FibroblastO,features = c("ELMOD2"), split.by = "tissue", ncol= 3, order = TRUE)

DotPlot(FibroblastO, features = c("NGFR"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
FeaturePlot(FibroblastO,features = c("NGFR","PAX9"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral,features = c("CCL5","ELMOD2"), split.by = "celltype.tissue", ncol= 3, order = TRUE)

DotPlot(FibroblastO, features = c("TGFB1","LOXL2","CCL2","CTGF","POSTN","THBS2"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(FibroblastO, features = c("TGFB1"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

DotPlot(oral, features = c("CTGF","POSTN","LOXL2","THBS2","ELMOD2"), group.by = "celltype.tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(oral, features = c("CXCL9","ELMOD2"), group.by = "celltype.tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

DotPlot(oral, features = c("IL6","IL1A","CXCL9","CXCL10"), group.by = "celltype.tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

FeaturePlot(oral,features = c("CD52","CCL3","CXCL9"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral,features = c("CCL5","ELMOD2"), split.by = "celltype.tissue", ncol= 3, order = TRUE)

DimPlot(oral, group.by = "cell.type")

FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10"), split.by = "tissue", ncol= 7, order = TRUE)

FibroblastO <- RunPCA(FibroblastO, verbose = TRUE)
ElbowPlot(FibroblastO, ndims = 50)
FibroblastO <- RunUMAP(FibroblastO, dims = 1:50, verbose = FALSE)

FibroblastO <- FindNeighbors(FibroblastO, dims = 1:50)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{FibroblastO <- FindClusters(FibroblastO, resolution = resolution)
}

library(clustree)

pdf('FibroblastO_clustree_seurat.pdf', width = 30, height = 20)
clustree(FibroblastO,node_size_range=c(10,20), node_text_size = 8)
dev.off()

Idents(object=FibroblastO) <- "SCT_snn_res.0.1"

DimPlot(FibroblastO, group.by = "SCT_snn_res.0.3")

FeaturePlot(all,features = c("HMOX1"), split.by = "tissue", ncol= 7, order = TRUE, pt.size = 0.1)

FeaturePlot(all,features = c("NT5E","THY1","ENG"), split.by = "tissue", ncol= 7, order = TRUE)
FeaturePlot(all,features = c("CD34","PTPRC"), split.by = "tissue", ncol= 7, order = TRUE)
FeaturePlot(all,features = c("FDCSP","TPSAB1","F13A1"), split.by = "tissue", ncol= 7, order = TRUE)

FeaturePlot(all,features = c("CXCL13","CTHRC1","POSTN"), split.by = "tissue", ncol= 7, order = TRUE, pt.size = 0.1)

FeaturePlot(all,features = c("CDKN2A"), split.by = "tissue", ncol= 7, order = TRUE, pt.size = 0.1)

FeaturePlot(all,features = c("CXCL13"), split.by = "tissue", ncol= 7, pt.size = 0.1, order = FALSE)
FeaturePlot(all,features = c("CTHRC1"), split.by = "tissue", ncol= 7, pt.size = 0.1, order = FALSE)
FeaturePlot(all,features = c("POSTN"), split.by = "tissue", ncol= 7, pt.size = 0.1, order = FALSE)

FeaturePlot(all,features = c("FDCSP"), split.by = "tissue", ncol= 7, pt.size = 0.1, order = FALSE)
FeaturePlot(all,features = c("TPSAB1"), split.by = "tissue", ncol= 7, pt.size = 0.1, order = FALSE)
FeaturePlot(all,features = c("F13A1"), split.by = "tissue", ncol= 7, pt.size = 0.1, order = FALSE)

FeaturePlot(all,features = c("OCT4","SSEA4","STRO1","SOX2","TP64","TP75"), split.by = "tissue", ncol= 7, order = TRUE)

pdf('PERICYTE.pdf', width = 19, height = 5)
FeaturePlot(all,features = c("PDGFRB","MCAM","CSPG4","KCNJ8","ABCC9","MYH11","HIGD1B","KCNA5","ACTA2","CSPG4","PLN","RERGL","CD248"), ncol= 6, order = TRUE)
dev.off()

FeaturePlot(allP,features = c("PDGFRB","MCAM","CSPG4","KCNJ8","ABCC9","MYH11","HIGD1B","KCNA5"), ncol= 5, order = TRUE)
#OTHER PERICYTE MARKERS: "KCNK3",

FeaturePlot(FibroblastO,features = c("KRT1","KRT2","KRT13","SLPI","KRT17","ANXA1","MYL6","ERP29"), split.by = "tissue", ncol= 7, order = TRUE)

FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), split.by = "tissue", ncol= 7, order = TRUE)
FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), split.by = "tissue", ncol= 7, order = FALSE)

FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), ncol= 7, order = TRUE)
FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), ncol= 7, order = TRUE)


FeaturePlot(FibroblastO,features = c("CXCL9","CXCL10","PITX1","WNT5A"), ncol= 2, order = TRUE)
FeaturePlot(FibroblastO,features = c("WNT5A"), ncol= 3, order = TRUE,split.by = "tissue")
FeaturePlot(FibroblastO,features = c("WNT5A"), ncol= 3, order = FALSE,split.by = "tissue")
FeaturePlot(FibroblastO,features = c("PITX1"), ncol= 3, order = TRUE,split.by = "tissue")
FeaturePlot(FibroblastO,features = c("PITX1"), ncol= 3, order = FALSE,split.by = "tissue")

Idents(object=FibroblastO) <- "SCT_snn_res.0.3"

saveRDS(FibroblastO, file = "FibroblastO.RDS")

FibroblastO.markers_0.3_res <- FindAllMarkers(FibroblastO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(openxlsx)
write.xlsx(FibroblastO.markers_0.3_res, file = 'FibroblastO.markers_0.3_res.xlsx', rowNames=TRUE)


FibroblastO.markers_0.3_res %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pdf('FibroblastO_markers_heatmap_top10.pdf', width = 8, height = 10)
DoHeatmap(FibroblastO, features = top10$gene) + NoLegend()
dev.off()

VlnPlot(FibroblastO, features = c("CXCL9","CXCL10","PITX1","WNT5A"), group.by = "SCT_snn_res.0.3", split.by = "tissue")
VlnPlot(FibroblastO, features = c("CXCL10"), group.by = "SCT_snn_res.0.3", split.by = "tissue")

FibroblastO$CXCL <- "no"

cxcl9pos <- WhichCells(FibroblastO, expression = CXCL9 > 0.5)
cxcl10pos <- WhichCells(FibroblastO, expression = CXCL10 > 0.5)


Idents(FibroblastO) <- "CXCL"
Idents(FibroblastO)
FibroblastO <- SetIdent(FibroblastO, cells =c(cxcl9pos,cxcl10pos), value = 'yes')
Idents(FibroblastO, cells =cxcl10pos) <- "yes"
FibroblastO$CXCL <- Idents(FibroblastO)
table(Idents(FibroblastO))


FibroblastO.CXCL9.10.markers <- FindMarkers(FibroblastO, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(FibroblastO.CXCL9.10.markers, file = 'FibroblastO.CXCL9.10.DEG.xlsx', rowNames=TRUE)



Idents(all)
table(all$cell.type)

all$celltype.tissue <- paste(Idents(all), all$tissue, sep = "_")
all$celltype <- Idents(all)
Idents(all) <- "celltype.tissue"


table(Idents(all))

count <-0
levels(Idents(all))
table(all$celltype) ###sanity check
for (ident in levels(Idents(all))) {
  if (length( WhichCells(all, ident = ident))==1) {
    tempmatrix<-as.matrix(GetAssayData(all, slot = "counts")[, WhichCells(all, ident = ident)])
  } else {
    tempmatrix<-as.matrix(Matrix::rowSums(GetAssayData(all, slot = "counts")[, WhichCells(all, ident = ident)]))
  }
  count <- count + 1
  if (match(ident,levels(Idents(all)))==1) {
    finalmatrix<-tempmatrix
  } else {
    finalmatrix<-cbind(finalmatrix,tempmatrix)
  }
}
colnames(finalmatrix)<-as.character(levels(Idents(all)))
write.table(finalmatrix, file='~/Mucosa/Dataset1/All-oral_Raw_Gene_Counts_per_Celltype-tissue.tsv', quote=TRUE, sep='\t', col.names = TRUE)
val1=as.numeric(5) ## MINCOUNTS
val2=as.numeric(1) ## MINSAMPLES
filter <- apply(finalmatrix, 1, function(x) length(x[x>val1])>=val2)
fmtxfiltered=finalmatrix[filter,]
write.table(as.data.frame(fmtxfiltered),file='~/Mucosa/Dataset1/All-oral_Raw_Gene_Counts_per_Celltype-tissue_filtered.tsv', quote=TRUE, sep='\t', col.names = TRUE)


Idents(all) <- "celltype"

table(Idents(all))

count <-0
levels(Idents(all))
table(all$celltype) ###sanity check
for (ident in levels(Idents(all))) {
  if (length( WhichCells(all, ident = ident))==1) {
    tempmatrix<-as.matrix(GetAssayData(all, slot = "counts")[, WhichCells(all, ident = ident)])
  } else {
    tempmatrix<-as.matrix(Matrix::rowSums(GetAssayData(all, slot = "counts")[, WhichCells(all, ident = ident)]))
  }
  count <- count + 1
  if (match(ident,levels(Idents(all)))==1) {
    finalmatrix<-tempmatrix
  } else {
    finalmatrix<-cbind(finalmatrix,tempmatrix)
  }
}
colnames(finalmatrix)<-as.character(levels(Idents(all)))
write.table(finalmatrix, file='~/Mucosa/Dataset1/All-oral_Raw_Gene_Counts_per_Celltype.tsv', quote=TRUE, sep='\t', col.names = TRUE)
val1=as.numeric(5) ## MINCOUNTS
val2=as.numeric(1) ## MINSAMPLES
filter <- apply(finalmatrix, 1, function(x) length(x[x>val1])>=val2)
fmtxfiltered=finalmatrix[filter,]
write.table(as.data.frame(fmtxfiltered),file='~/Mucosa/Dataset1/All-oral_Raw_Gene_Counts_per_Celltype_filtered.tsv', quote=TRUE, sep='\t', col.names = TRUE)


Idents(all) <- "tissue"
table(Idents(all))

oral<-subset(all, idents = c("Buccal","Gingiva","Perio"))

Idents(oral) <- "celltype"

table(Idents(oral))

count <-0
levels(Idents(oral))
table(oral$celltype) ###sanity check
for (ident in levels(Idents(oral))) {
  if (length( WhichCells(oral, ident = ident))==1) {
    tempmatrix<-as.matrix(GetAssayData(oral, slot = "counts")[, WhichCells(oral, ident = ident)])
  } else {
    tempmatrix<-as.matrix(Matrix::rowSums(GetAssayData(oral, slot = "counts")[, WhichCells(oral, ident = ident)]))
  }
  count <- count + 1
  if (match(ident,levels(Idents(oral)))==1) {
    finalmatrix<-tempmatrix
  } else {
    finalmatrix<-cbind(finalmatrix,tempmatrix)
  }
}
colnames(finalmatrix)<-as.character(levels(Idents(oral)))
write.table(finalmatrix, file='~/Mucosa/Dataset1/oral_Raw_Gene_Counts_per_Celltype.tsv', quote=TRUE, sep='\t', col.names = TRUE)
val1=as.numeric(5) ## MINCOUNTS
val2=as.numeric(1) ## MINSAMPLES
filter <- apply(finalmatrix, 1, function(x) length(x[x>val1])>=val2)
fmtxfiltered=finalmatrix[filter,]
write.table(as.data.frame(fmtxfiltered),file='~/Mucosa/Dataset1/oral_Raw_Gene_Counts_per_Celltype_filtered.tsv', quote=TRUE, sep='\t', col.names = TRUE)




saveRDS(oral, file = "oral.RDS")
library(future)
plan("multicore", workers = 1)
options(future.globals.maxSize = 32000 * 1024^2)

oral <- SCTransform(oral, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
gc()
oral <- RunPCA(oral, verbose = TRUE)
ElbowPlot(oral, ndims = 50)
oral <- RunUMAP(oral, dims = 1:50, verbose = FALSE)

oral <- FindNeighbors(oral, dims = 1:50)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{oral <- FindClusters(oral, resolution = resolution)
}

library(clustree)

pdf('oral_clustree_seurat.pdf', width = 30, height = 20)
clustree(oral,node_size_range=c(10,20), node_text_size = 8)
dev.off()

oral$SCT_snn_res.0.1

Idents(object=oral) <- "SCT_snn_res.0.1"
Idents(oral)

saveRDS(oral, file = "oral_normalized_and_clustered_umap.RDS")

FeaturePlot(oral,features = c("PDGFRB","MCAM","CSPG4","KCNJ8","ABCC9","MYH11","HIGD1B","KCNA5"), split.by = "tissue", ncol= 3, order = TRUE)


FeaturePlot(oral,features = c("GPX1","GPX2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral,features = c("GPX3","GPX4"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral,features = c("HMOX1","NFE2L2"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(all,features = c("TNMD","FBN2","TBX3","PRKCDBP","LHFP","PTRF","GPX1","SEPP1","ATP5L","TCEB2","USMG5","ALDOA","NGFRAP1","SHFM1","IGFBP2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("GPX1","SEPP1","ATP5L"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("TCEB2","USMG5","ALDOA"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("NGFRAP1","SHFM1","TBX3"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(all,features = c("PRKCDBP","LHFP"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("STMN2","PRKCDBP"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("TCEB1","ELOA"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("ELOA","ELOB","ELOC"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("IGF2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("STMN2"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(all,features = c("CLEC3B","HMOX1","FAM180B"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("FOXQ1","IRX2","COL4A4"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("COL4A1","COL4A2","COL4A3","COL4A4"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("PITX1","WNT5A"), split.by = "tissue", ncol= 3, order = TRUE)
DimPlot(all, group.by = "tissue")
FeaturePlot(all,features = c("PRKCDBP","GPX1"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("IGF2","WNT5A"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("STMN2","ELOA"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("STMN2","FOXQ1","ELOA"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(all,features = c("PRKCDBP","GPX1","IGF2","WNT5A"), split.by = "tissue", ncol= 3, order = TRUE)


DimPlot(oral, group.by = "celltype")
DimPlot(oral, group.by = "celltype")
DimPlot(oral, group.by = "tissue")

DimPlot(all, group.by = "celltype")
DimPlot(all, group.by = "tissue")



oral<- RenameIdents(oral, 'myofibroblast'='pericyte')

table(Idents(oral))

oral$celltype <-Idents(oral)

Idents(oral) <- "celltype"
table(Idents(oral))
oral_epi<-subset(oral, idents = c("epithelium"))

oral_epi <- SCTransform(oral_epi, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
gc()
oral_epi <- RunPCA(oral_epi, verbose = TRUE)
ElbowPlot(oral_epi, ndims = 50)
oral_epi <- RunUMAP(oral_epi, dims = 1:50, verbose = FALSE)

oral_epi <- FindNeighbors(oral_epi, dims = 1:50)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{oral_epi <- FindClusters(oral_epi, resolution = resolution)
}

library(clustree)

pdf('oral_epi_clustree_seurat.pdf', width = 30, height = 20)
clustree(oral_epi,node_size_range=c(10,20), node_text_size = 8)
dev.off()

oral_epi$SCT_snn_res.0.1

Idents(object=oral_epi) <- "SCT_snn_res.0.1"
Idents(oral_epi)


saveRDS(oral_epi, file = "oral_epi_normalized_and_clustered_umap.RDS")


FeaturePlot(oral_epi,features = c("KRT1","SLPI","KRT17","ANXA1","MYL6","ERP29"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("KRT1","KRT2","KRT13","SLPI","KRT17","ANXA1","MYL6","ERP29"), split.by = "celtype", ncol= 7, order = TRUE)

oral_epi<-subset(oral, idents = c("epithelium"))

oral_epi <- SCTransform(oral_epi, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
gc()
oral_epi <- RunPCA(oral_epi, verbose = TRUE)
ElbowPlot(oral_epi, ndims = 50)
oral_epi <- RunUMAP(oral_epi, dims = 1:7, verbose = FALSE)

oral_epi <- FindNeighbors(oral_epi, dims = 1:7)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{oral_epi <- FindClusters(oral_epi, resolution = resolution)
}

library(clustree)

pdf('oral_epi_clustree_seurat.pdf', width = 30, height = 20)
clustree(oral_epi,node_size_range=c(10,20), node_text_size = 8)
dev.off()

oral_epi$SCT_snn_res.0.1

Idents(object=oral_epi) <- "SCT_snn_res.0.1"
Idents(oral_epi)

FeaturePlot(oral_epi,features = c("KRT1","SLPI","KRT17","ANXA1","MYL6","ERP29"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("ODAM","SLPI","KRT17"), split.by = "tissue", ncol= 3, order = TRUE)
DimPlot(oral_epi, group.by = "tissue")
DimPlot(oral_epi)

FeaturePlot(oral_epi,features = c("BMP2"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(oral_epi,features = c("BMP2","MMP7"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("MMP13","SLC7A11"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(oral_epi,features = c("GPX1","GPX2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("GPX3","GPX4"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("HMOX1","NFE2L2"), split.by = "tissue", ncol= 3, order = TRUE)

FeaturePlot(oral_epi,features = c("CSF3","LAMC2"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
                                  "BMP2","MMP7","MMP13","SLC7A11"), split.by = "tissue", ncol= 3, order = TRUE)


library(ggplot2)

pdf('JE_DotPlot.pdf', width = 30, height = 20)
#DotPlot(oral_epi, features = c("DST","HMGN3","GLTSCR2","KRT17","SLPI","ODAM","SPRR2F","S100A7","IL36G","LCN2","DNASE1L3","CST6","PRSS22","ERO1A","CSF3","ODC1","HBEGF","TNFAIP3"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(oral_epi, features = c("DST","HMGN3","GLTSCR2","ALDH3A1","KRT17","ODAM","SPRR2F","S100A7","IL36G","LCN2","MMP12","C4orf26","SAA2","CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
                               "BMP2","MMP7","MMP13","SLC7A11"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
#DotPlot(oral_epi, features = c("DST","HMGN3","GLTSCR2","ALDH3A1","KRT17","ODAM","SLPI","SPRR2F","S100A7","IL36G","LCN2","MMP12","C4orf26","SAA2","CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
"BMP2","MMP7","MMP13","SLC7A11"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
#DotPlot(oral_epi, features = c("DST","HMGN3","GLTSCR2","ALDH3A1","KRT17","ODAM","SLPI","SPRR2F","S100A7","IL36G","LCN2","MMP12","C4orf26","SAA2","CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","QPCT","VMO1","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
#                               "MMP7","BMP2","MMP13","SLC7A11"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
#DotPlot(oral_epi, features = c("KRT15","DST","HMGN3","GLTSCR2","KRT17","SLPI","ODAM","SPRR2F","S100A7","IL36G","LCN2","DNASE1L3","CST6","PRSS22","ERO1A","CSF3","ODC1","UBE2S","TNFAIP3"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

DotPlot(oral_epi, features = c("KRT17","ODAM","SPRR2F","S100A7","IL36G","LCN2","MMP12","C4orf26","SAA2","CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
                               "BMP2","MMP7","MMP13","SLC7A11","DST","HMGN3","GLTSCR2","ALDH3A1"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

DotPlot(oral_epi, features = c("ITGA1","ITGA2","ITGA3","ITGA4","ITGA5","ITGA6","ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
FeaturePlot(oral_epi,features = c("COL17A1","DST","PLEC","ITGA6"), split.by = "tissue", ncol= 3, order = TRUE)

pdf('JE_umap3.pdf', width = 4.4, height = 3.5)
DimPlot(oral_epi, group.by = "tissue", shuffle = "TRUE")
dev.off()

FeaturePlot(oral_epi,features = c("ODAM"), ncol= 3, order = TRUE)

#UP IN PERIO EPI
FeaturePlot(oral_epi,features = c("PLAUR","TNFRSF21"), split.by = "tissue", ncol= 3, order = TRUE)
#UP IN HEALTHY EPI
FeaturePlot(oral_epi,features = c("CXCL8","CD55"), split.by = "tissue", ncol= 3, order = TRUE)
FeaturePlot(oral_epi,features = c("CLDN1","ST14"), split.by = "tissue", ncol= 3, order = TRUE)

saveRDS(oral_epi, file = "oral_epi_normalized_and_clustered_umap.RDS")



oral_epi$krt17 <- "no"

oral_epi$krt17

oral_epi

krt171pos <- WhichCells(oral_epi, expression = KRT17 > 0.5) ## 1364 cells

oral_epi ##total number of cells 7873


Idents(oral_epi) <- "krt17"

Idents(oral_epi)

oral_epi <- SetIdent(oral_epi, cells =krt171pos, value = 'yes')

table(Idents(oral_epi), oral_epi$tissue)

#KRT17
#    Buccal Gingiva Perio
#yes     52     875   437
#no    1722    4326   461

oral_epi$krt17 <- Idents(oral_epi)

table(Idents(oral_epi))
library(openxlsx)

oral_epi.krt17.markers <- FindMarkers(oral_epi, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")

write.xlsx(oral_epi.krt17.markers, file = 'oral_epi.krt17.DEG.xlsx', rowNames=TRUE)

oral_epi$tissue.krt17pos <- paste(oral_epi$tissue, Idents(oral_epi), sep = "_")
Idents(oral_epi) <- "tissue.krt17pos"
Idents(oral_epi)

oralkrt17pos_gingiva_vs_perio.markers <- FindMarkers(oral_epi, ident.1 = "Gingiva_yes", ident.2 = "Perio_yes", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(oralkrt17pos_gingiva_vs_perio.markers, file = 'oralkrt17pos_gingiva_vs_perio.DEG.xlsx', rowNames=TRUE)

oralNOkrt17pos_gingiva_vs_perio.markers <- FindMarkers(oral_epi, ident.1 = "Gingiva_no", ident.2 = "Perio_no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(oralNOkrt17pos_gingiva_vs_perio.markers, file = 'oralNOkrt17pos_gingiva_vs_perio.DEG.xlsx', rowNames=TRUE)


#####ODAM
#KRT17
#    Buccal Gingiva Perio
#yes     52     875   437
#no    1722    4326   461

oral_epi$odam <- "no"

oral_epi$odam

oral_epi

odam1pos <- WhichCells(oral_epi, expression = ODAM > 0.5) ## 1364 cells

oral_epi ##total number of cells 7873


Idents(oral_epi) <- "odam"

Idents(oral_epi)

oral_epi <- SetIdent(oral_epi, cells =odam1pos, value = 'yes')

table(Idents(oral_epi), oral_epi$tissue)

#KRT17
#    Buccal Gingiva Perio
#yes     52     875   437
#no    1722    4326   461

#ODAM        Buccal Gingiva Perio
#Buccal_no     1774       0     0
#Gingiva_no       0    4280     0
#Gingiva_yes      0     921     0
#Perio_no         0       0   625
#Perio_yes        0       0   273

oral_epi$odam <- Idents(oral_epi)

table(Idents(oral_epi))
library(openxlsx)

oral_epi.odam.markers <- FindMarkers(oral_epi, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")

write.xlsx(oral_epi.odam.markers, file = 'oral_epi.odam.DEG.xlsx', rowNames=TRUE)

oral_epi$tissue.odampos <- paste(oral_epi$tissue, Idents(oral_epi), sep = "_")
Idents(oral_epi) <- "tissue.odampos"
Idents(oral_epi)

oralodampos_gingiva_vs_perio.markers <- FindMarkers(oral_epi, ident.1 = "Gingiva_yes", ident.2 = "Perio_yes", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(oralodampos_gingiva_vs_perio.markers, file = 'oralodampos_gingiva_vs_perio.DEG.xlsx', rowNames=TRUE)

oralNOodampos_gingiva_vs_perio.markers <- FindMarkers(oral_epi, ident.1 = "Gingiva_no", ident.2 = "Perio_no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
write.xlsx(oralNOodampos_gingiva_vs_perio.markers, file = 'oralNOodampos_gingiva_vs_perio.DEG.xlsx', rowNames=TRUE)


### ISOLATE ODAM+ JE CELLS??

Idents(oral_epi)
offendercell <- WhichCells(oral_epi, expression = ODAM > 0.5)
JE <- subset(oral_epi, cells = offendercell, invert = FALSE)

FeaturePlot(JE,features = c("ODAM"), split.by = "tissue", order = TRUE)


JE <- SCTransform(JE, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
gc()
JE <- RunPCA(JE, verbose = TRUE)
ElbowPlot(JE, ndims = 30)
JE <- RunUMAP(JE, dims = 1:10, verbose = FALSE)

JE <- FindNeighbors(JE, dims = 1:10)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
{JE <- FindClusters(JE, resolution = resolution)
}

DimPlot(JE, group.by = "tissue")
DimPlot(JE, group.by = "SCT_snn_res.0.01")
DimPlot(JE)

FeaturePlot(JE,features = c("ODAM"), split.by = "tissue", order = TRUE)

pdf('JEonly_DotPlot.pdf', width = 30, height = 20)
#DotPlot(oral_epi, features = c("DST","HMGN3","GLTSCR2","KRT17","SLPI","ODAM","SPRR2F","S100A7","IL36G","LCN2","DNASE1L3","CST6","PRSS22","ERO1A","CSF3","ODC1","HBEGF","TNFAIP3"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(JE, features = c("KRT17","ODAM","SPRR2F","S100A7","IL36G","LCN2","MMP12","C4orf26","SAA2","CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
                         "BMP2","MMP7","MMP13","SLC7A11"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
#DotPlot(JE, features = c("DST","HMGN3","GLTSCR2","ALDH3A1","KRT17","ODAM","SPRR2F","S100A7","IL36G","LCN2","MMP12","C4orf26","SAA2","CD74","COL3A1","COL1A2","COL1A1","DNASE1L3","CST6","PRSS22","ERO1A","TFRC","CSF3","ODC1","HBEGF","TNFAIP3","LAMC2",
#                               "BMP2","MMP7","MMP13","SLC7A11"), group.by = "tissue") + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

DotPlot(JE, features = c("ODAM"), group.by = "tissue", scale.min = 0, scale.max = 99, col.min = 0, col.max = 1) + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))


JE$SCT_snn_res.0.1

Idents(object=JE) <- "SCT_snn_res.0.1"
Idents(JE)

saveRDS(JE, file = "JE_normalized_and_clustered_umap.RDS")

##########MSC

#MSC STEM CELL MARKERS - ALSO ,"PAX9"
FeaturePlot(FibroblastO,features = c("NT5E","ENG","THY1","CD44"), split.by = "tissue", ncol= 3, order = TRUE)

#NEGATIVE MSC MARKER: CD14, CD19, CD34, CD45 (PTRC)
FeaturePlot(FibroblastO,features = c("CD14","CD19","CD34","PTPRC"), split.by = "tissue", ncol= 3, order = TRUE)

FibroblastO$MSC <- "no"
#MSC positive markers
nt5epos <- WhichCells(FibroblastO, expression = NT5E > 0.5)
engpos <- WhichCells(FibroblastO, expression = ENG > 0.5)
thy1pos <- WhichCells(FibroblastO, expression = THY1 > 0.5)
cd44pos <- WhichCells(FibroblastO, expression = CD44 > 0.5)

#MSC negative markers
cd14pos <- WhichCells(FibroblastO, expression = CD14 > 0.5)
cd19pos <- WhichCells(FibroblastO, expression = CD19 > 0.5)
cd34pos <- WhichCells(FibroblastO, expression = CD34 > 0.5)
ptprcpos <- WhichCells(FibroblastO, expression = PTPRC > 0.5)

Idents(FibroblastO) <- "MSC"

Idents(FibroblastO)


FibroblastO <- SetIdent(FibroblastO, cells =c(nt5epos,engpos,thy1pos,cd44pos), value = 'yes')
FibroblastO <- SetIdent(FibroblastO, cells =c(cd14pos,cd19pos,cd34pos,ptprcpos), value = 'no')

Idents(FibroblastO, cells =c(nt5epos,engpos,thy1pos,cd44pos)) <- "yes"

FibroblastO$MSC <- Idents(FibroblastO)

table(Idents(FibroblastO))


DefaultAssay(FibroblastO)<-"RNA"

FibroblastO$MSC.tissue <- paste(Idents(FibroblastO), FibroblastO$tissue, sep = "_")
FibroblastO$MSC <- Idents(FibroblastO)
Idents(FibroblastO) <- "MSC.tissue"
Idents(FibroblastO)
FibroblastO$MSC.tissue

fibro_DE_genes<-FindMarkers(FibroblastO, ident.1 = c("yes_Buccal","yes_Gingiva","yes_Perio"), ident.2 = c("no_Buccal","no_Gingiva","no_Perio"), only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.01, test.use="MAST")
write.xlsx(fibro_DE_genes, file = 'MSC_DEG.xlsx', rowNames=TRUE)

fibro_DE_genes2<-FindMarkers(FibroblastO, ident.1 = c("yes_Buccal","yes_Gingiva"), ident.2 = c("yes_Gingiva","yes_Perio"), only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.01, test.use="MAST")
write.xlsx(fibro_DE_genes2, file = 'MSC_DEG.xlsx', rowNames=TRUE)

FibroblastO.MSC.markers <- FindMarkers(FibroblastO, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")
library(openxlsx)
write.xlsx(FibroblastO.MSC.markers, file = 'FibroblastO.MSC.DEG.xlsx', rowNames=TRUE)

###version 2
Idents(FibroblastO, cells =c(nt5epos,engpos,thy1pos,cd44pos,cd14pos,cd19pos,cd34pos,ptprcpos)) <- "yes"

FibroblastO$MSC <- Idents(FibroblastO)

table(Idents(FibroblastO))


FibroblastO.MSC.markers <- FindMarkers(FibroblastO, ident.1 = "yes", ident.2 = "no", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")

write.xlsx(FibroblastO.MSC.markers, file = 'FibroblastO.MSCwithNegative.DEG.xlsx', rowNames=TRUE)
