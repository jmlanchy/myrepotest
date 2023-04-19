library(Seurat)
library(tidyverse)

filepath <- "C:/Users/jmlan/Documents/Projects_Mark_Grimes/single_cell_nuclei_RNAseq/mg_snrna_march_2023/p1d87/outs/"
sobj87.data <- Read10X_h5(paste0(filepath, "filtered_feature_bc_matrix.h5", sep = ""))
sobj87 <- CreateSeuratObject(sobj87.data[['Gene Expression']])
sobj87[["HTO"]] <- CreateAssayObject(sobj87.data[["Antibody Capture"]])

sobj87 <- NormalizeData(sobj87)                                
sobj87 <- FindVariableFeatures(sobj87)                          
sobj87 <- ScaleData(sobj87, features = VariableFeatures(sobj87))
sobj87 <- RunPCA(sobj87, features = VariableFeatures(sobj87))

RunTSNE(
  sobj87,
  dims = 1:20,
  perplexity = 30,
  reduction = 'pca') -> sobj87

FindNeighbors(sobj87, dims=1:20, reduction = "pca") -> sobj87
FindClusters(sobj87, resolution = 0.1) -> sobj87

RunUMAP(
  sobj87,
  dims=1:20,
  reduction = 'pca') -> sobj87

sobj87.markers <- FindAllMarkers(sobj87, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj87.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_87
dev.new()
DoHeatmap(sobj87, features = top10_87$gene) + NoLegend()

# From Mark's
# Compare with ROC
# Alternative: the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect)
sobj87.markers.roc <- FindAllMarkers(sobj87, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, use = "roc")
sobj87.markers.roc %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_87roc
dev.new()
DoHeatmap(sobj87, features = top10_87roc$gene) + NoLegend()

# I run day1

filepath1 <- "C:/Users/jmlan/Documents/Projects_Mark_Grimes/single_cell_nuclei_RNAseq/mg_snrna_march_2023/p1d1/outs/"
sobj1.data <- Read10X_h5(paste0(filepath1, "filtered_feature_bc_matrix.h5", sep = ""))
sobj1 <- CreateSeuratObject(sobj1.data[['Gene Expression']])
sobj1[["HTO"]] <- CreateAssayObject(sobj1.data[["Antibody Capture"]])

sobj1 <- NormalizeData(sobj1)                                
sobj1 <- FindVariableFeatures(sobj1)                          
sobj1 <- ScaleData(sobj1, features = VariableFeatures(sobj1))
sobj1 <- RunPCA(sobj1, features = VariableFeatures(sobj1))

RunTSNE(
  sobj1,
  dims=1:20,
  perplexity=30,
  reduction = 'pca'
) -> sobj1

FindNeighbors(sobj1, dims=1:20, reduction = "pca") -> sobj1
FindClusters(sobj1,resolution = 0.1) -> sobj1

RunUMAP(
  sobj1,
  dims=1:20,
  reduction = 'pca'
) -> sobj1

UMAPPlot(sobj1)


sobj1.markers <- FindAllMarkers(sobj1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj1.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_1
dev.new()
DoHeatmap(sobj1, features = top10_1$gene) + NoLegend()
# Compare with ROC
# Alternative: the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect)
sobj1.markers.roc <- FindAllMarkers(sobj1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, use = "roc")
sobj1.markers.roc %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_1roc
dev.new()
DoHeatmap(sobj1, features = top10_1roc$gene) + NoLegend()

####### List of markers ################
### I first had a look at markers with a simple pattern, like ^COL.
### Looking at the list obtained thru the simple pattern screening,
### I wrote a regular expression filter to exclude the markers having apparently no links with
### the expected markers' characteristics.

# ribosomal protein genes
PercentageFeatureSet(sobj1, pattern = "^RP[LS]") -> sobj1$percent.ribosomal
PercentageFeatureSet(sobj87, pattern = "^RP[LS]") -> sobj87$percent.ribosomal
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.ribosomal"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.ribosomal"))

# collagen protein genes
PercentageFeatureSet(sobj1, pattern = "^COL[0-9]{1,2}A[0-9]$") -> sobj1$percent.collagen
PercentageFeatureSet(sobj87, pattern = "^COL[0-9]{1,2}A[0-9]$") -> sobj87$percent.collagen
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.collagen"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.collagen"))
collagen.markers <- grep("^COL[0-9]{1,2}A[0-9]$", rownames(sobj1@assays$RNA@counts), value = TRUE)
VlnPlot(sobj1, features = collagen.markers[1:16])
VlnPlot(sobj87, features =  collagen.markers[1:16])

# WNT protein genes
PercentageFeatureSet(sobj1, pattern = "^WNT[0-9]{1,2}$|^WNT[0-9]{1,2}[AB]$") -> sobj1$percent.wnt
PercentageFeatureSet(sobj87, pattern = "^WNT[0-9]{1,2}$|^WNT[0-9]{1,2}[AB]$") -> sobj87$percent.wnt
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.wnt"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.wnt"))
# from now on, I integrated the grep step in the VlnPlot command
VlnPlot(sobj1, features =  grep("^WNT[0-9]{1,2}$|^WNT[0-9]{1,2}[AB]$", rownames(sobj1@assays$RNA@counts), value = TRUE))
VlnPlot(sobj87, features =  grep("^WNT[0-9]{1,2}$|^WNT[0-9]{1,2}[AB]$", rownames(sobj1@assays$RNA@counts), value = TRUE))

# PDGF protein genes
PercentageFeatureSet(sobj1, pattern = "^PDGF") -> sobj1$percent.pdgf
PercentageFeatureSet(sobj87, pattern = "^PDGF") -> sobj87$percent.pdgf
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.pdgf"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.pdgf"))
VlnPlot(sobj1, features =  grep("^PDGF", rownames(sobj1@assays$RNA@counts), value = TRUE))
VlnPlot(sobj87, features =  grep("^PDGF", rownames(sobj1@assays$RNA@counts), value = TRUE))

# NOTCH protein genes
PercentageFeatureSet(sobj1, pattern = "^NOTCH") -> sobj1$percent.notch
PercentageFeatureSet(sobj87, pattern = "^NOTCH") -> sobj87$percent.notch
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.notch"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.notch"))
VlnPlot(sobj1, features =  grep("^NOTCH", rownames(sobj1@assays$RNA@counts), value = TRUE))
VlnPlot(sobj87, features =  grep("^NOTCH", rownames(sobj1@assays$RNA@counts), value = TRUE))

# SOX protein genes
PercentageFeatureSet(sobj1, pattern = "^SOX[0-9]{1,2}$") -> sobj1$percent.sox
PercentageFeatureSet(sobj87, pattern = "^SOX[0-9]{1,2}$") -> sobj87$percent.sox
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.sox"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.sox"))
VlnPlot(sobj1, features =  grep("^SOX[0-9]{1,2}$",rownames(sobj1@assays$RNA@counts), value = TRUE))
VlnPlot(sobj87, features =  grep("^SOX[0-9]{1,2}$",rownames(sobj1@assays$RNA@counts), value = TRUE))

# BMP protein genes
PercentageFeatureSet(sobj1, pattern = "^BMP[0-9]{1,2}$|^BMP8[AB]$|^BMPER$|^BMPR1[AB]$|^BMPR2$") -> sobj1$percent.sox
PercentageFeatureSet(sobj87, pattern = "^BMP[0-9]{1,2}$|^BMP8[AB]$|^BMPER$|^BMPR1[AB]$|^BMPR2$") -> sobj87$percent.sox
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.sox"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.sox"))
VlnPlot(sobj1, features =  grep("^BMP[0-9]{1,2}$|^BMP8[AB]$|^BMPER$|^BMPR1[AB]$|^BMPR2$",rownames(sobj1@assays$RNA@counts), value = TRUE))
VlnPlot(sobj87, features =  grep("^BMP[0-9]{1,2}$|^BMP8[AB]$|^BMPER$|^BMPR1[AB]$|^BMPR2$",rownames(sobj1@assays$RNA@counts), value = TRUE))

# TGFB protein genes
PercentageFeatureSet(sobj1, pattern = "^TGFBI$|^TGFB[1-3]$|^TGFBR[1-3]$") -> sobj1$percent.tgfb
PercentageFeatureSet(sobj87, pattern = "^TGFBI$|^TGFB[1-3]$|^TGFBR[1-3]$") -> sobj87$percent.tgfb
VlnPlot(sobj1, features = c("nFeature_RNA", "nCount_RNA", "percent.tgfb"))
VlnPlot(sobj87, features = c("nFeature_RNA", "nCount_RNA", "percent.tgfb"))
VlnPlot(sobj1, features = grep("^TGFBI$|^TGFB[1-3]$|^TGFBR[1-3]$",rownames(sobj1@assays$RNA@counts), value = TRUE))
VlnPlot(sobj87, features = grep("^TGFBI$|^TGFB[1-3]$|^TGFBR[1-3]$",rownames(sobj1@assays$RNA@counts), value = TRUE))

sobj.combined <- merge(sobj1, y = sobj87, add.cell.ids = c("day01", "day87"), project = "cartilage")
sobj.combined
sobj.combined <- NormalizeData(sobj.combined)                                
sobj.combined <- FindVariableFeatures(sobj.combined)                          
sobj.combined <- ScaleData(sobj.combined, features = VariableFeatures(sobj.combined))
sobj.combined <- RunPCA(sobj.combined, features = VariableFeatures(sobj.combined))

RunTSNE(
  sobj.combined,
  dims=1:20,
  perplexity=30,
  reduction = 'pca'
) -> sobj.combined

FindNeighbors(sobj.combined, dims=1:20, reduction = "pca") -> sobj.combined
FindClusters(sobj.combined,resolution = 0.1) -> sobj.combined

RunUMAP(
  sobj.combined,
  dims=1:20,
  reduction = 'pca'
) -> sobj.combined

UMAPPlot(sobj.combined)
TSNEPlot(sobj.combined)

# here I load and process the hashed data

filepath2 <- "C:/Users/jmlan/Documents/Projects_Mark_Grimes/single_cell_nuclei_RNAseq/mg_snrna_march_2023/hashed_lane/outs/"
sobj.data <- Read10X_h5(paste0(filepath2, "filtered_feature_bc_matrix.h5", sep = ""))
sobj <- CreateSeuratObject(sobj.data[['Gene Expression']])
sobj[["HTO"]] <- CreateAssayObject(sobj.data[["Antibody Capture"]])

sobj <- NormalizeData(sobj)                                
sobj <- FindVariableFeatures(sobj)                          
sobj <- ScaleData(sobj, features = VariableFeatures(sobj))
sobj <- RunPCA(sobj, features = VariableFeatures(sobj))

RunTSNE(
  sobj,
  dims=1:20,
  perplexity=30,
  reduction = 'pca'
) -> sobj

FindNeighbors(sobj, dims=1:20, reduction = "pca") -> sobj
FindClusters(sobj,resolution = 0.1) -> sobj

RunUMAP(
  sobj,
  dims=1:20,
  reduction = 'pca'
) -> sobj

PercentageFeatureSet(sobj, pattern = "^COL[0-9]{1,2}A[0-9]$") -> sobj$percent.collagen
PercentageFeatureSet(sobj.combined, pattern = "^COL[0-9]{1,2}A[0-9]$") -> sobj.combined$percent.collagen
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.collagen"))
VlnPlot(sobj.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.collagen"))

PercentageFeatureSet(sobj, pattern = "^PDGF") -> sobj$percent.pdgf
PercentageFeatureSet(sobj.combined, pattern = "^PDGF") -> sobj.combined$percent.pdgf
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.pdgf"))
VlnPlot(sobj.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.pdgf"))
VlnPlot(sobj, features =  grep("^PDGF",rownames(sobj@assays$RNA@counts), value = TRUE))
VlnPlot(sobj.combined, features =  grep("^PDGF",rownames(sobj.combined@assays$RNA@counts), value = TRUE))
VlnPlot(sobj, features = grep("^TGFBI$|^TGFB[1-3]$|^TGFBR[1-3]$",rownames(sobj@assays$RNA@counts), value = TRUE))
VlnPlot(sobj.combined, features = grep("^TGFBI$|^TGFB[1-3]$|^TGFBR[1-3]$",rownames(sobj.combined@assays$RNA@counts), value = TRUE))

table(sobj@meta.data$seurat.clusters)
neuromarkers <- c("NRXN1", "NRXN2", "NLGN1", "MDGA2", "GPHN", "DLG4", "NTRK1", "NTRK2", "NTRK3")
VlnPlot(sobj, features = neuromarkers)
VlnPlot(sobj.combined, features = neuromarkers)

day1cluster2.markers <- FindMarkers(sobj1, ident.1 = 2, ident.2 = c(0, 1), min.pct = 0.25)
day1cluster1.markers <- FindMarkers(sobj1, ident.1 = 1, ident.2 = c(0, 2), min.pct = 0.25)
day1cluster0.markers <- FindMarkers(sobj1, ident.1 = 0, ident.2 = c(1, 2), min.pct = 0.25)

# let's plot the protein.nodedata
top100cartilagemarkersday87 <- c("COL1A1", "COL3A1", "COL1A2", "COL6A3", "POSTN", "PRRX1", "COL5A1", "COL6A1", "COL12A1", "DCN", "BICC1", "COL6A2", "SPARC", "LUM", "BNC2", "ZFHX3", "MGP", "VCAN", "DLC1", "ADAM12", "FBLN1", "HMGA2", "TMEM132C", "RBMS3", "TWIST2", "TGFBR3", "PRKG1", "H19", "CNTN4", "RERG", "COL5A2", "FN1", "OGN", "SGCD", "GLIS1", "SERPINH1", "PDGFRA", "LAMA4", "ARHGAP24", "DLK1", "SFRP2", "LAMB1", "CYP1B1", "PDZRN3", "DNM3OS", "EDNRA", "GREM1", "CDH11", "FLRT2", "NPR3", "PLPP3", "PMEPA1", "IGFBP4", "MFAP4", "SEMA3A", "ARHGAP28", "PDE7B", "TBX15", "COL26A1", "ST6GALNAC3", "ROR2", "NAV3", "PRICKLE1", "RRBP1", "STK32B", "MEOX2", "PDGFRB", "GXYLT2", "TWIST1", "EBF1", "COL16A1", "EEF1A1", "SEMA5A", "COLEC12", "DDR2", "SLIT3", "ADAMTS12", "EML4", "CDON", "MEST", "COL14A1", "ADAMTS9", "SMOC2", "SEC24D", "SVIL", "HMCN1", "COL5A3", "EGFL6", "MT-CO1", "RPS2", "KCNQ1OT1", "ROR1", "MT-CYB", "EFNA5", "RPL3", "RUNX1", "APCDD1", "LGALS1", "CCDC80")
ggplot(data = filter(geneid.protein.nodedata.rowz, (geneid %in% top100cartilagemarkersday87) == TRUE)) +
            geom_point(aes(x = NCSCmedians, y = orgmedians), size = 3, color = "blue") +
            geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
            geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
            geom_abline(slope = 1, linewidth = 1, color = "blue") +
            geom_label_repel(aes(x = NCSCmedians, y = orgmedians, label = geneid), size = 4) +
            theme_bw() +
            theme(axis.text.x = element_text(size = 24),
                  axis.text.y = element_text(size = 23),
                  axis.title.x = element_text(size = 24),
                  axis.title.y = element_text(size = 24)) +
            ggtitle('neuromarkers)') 
            # + xlim(-10,10)+ylim(-10,10)

ggplot(data = filter(geneid.protein.nodedata.rowz, (geneid %in% c("NGFR")) == TRUE)) +
          geom_point(aes(x = geneid, y = NCSCmedians), color = "blue", size = 5) +
          geom_point(aes(x = geneid, y = NCSCmax), color = "orange", size = 5) +
          ylab("NCSCmedians (blue) or NCSCmax (orange)")

# Let's isolate the top 500 markers of sobj1, sobj87 and sobj.combined (merged)
sobj1.markers.roc %>%
  group_by(cluster) %>%
  top_n(n = 500, wt = avg_log2FC) -> top500_1roc
sobj87.markers.roc %>%
  group_by(cluster) %>%
  top_n(n = 500, wt = avg_log2FC) -> top500_87roc
sobj.combined.markers.roc <- FindAllMarkers(sobj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, use = "roc")
sobj.combined.markers.roc %>%
  group_by(cluster) %>%
  top_n(n = 500, wt = avg_log2FC) -> top500_combinedroc

unique(top500_1roc$cluster)
[1] 0 1 2
Levels: 0 1 2
unique(top500_87roc$cluster)
[1] 0 1 2 3 4 5
Levels: 0 1 2 3 4 5
unique(top500_combinedroc$cluster)
[1] 0 1 2 3 4 5
Levels: 0 1 2 3 4 5

# let's test the intersection between the top500 genes of cluster0, day1 with cluster2, day1.
day1cluster0.list <- filter(top500_1roc, cluster == 0)$gene
day1cluster2.list <- filter(top500_1roc, cluster == 2)$gene

day87cluster0.list <- filter(top500_87roc, cluster == 0)$gene
day87cluster4.list <- filter(top500_87roc, cluster == 4)$gene

combined.cluster0.list <- filter(top500_combinedroc, cluster == 0)$gene
combined.cluster1.list <- filter(top500_combinedroc, cluster == 1)$gene
combined.cluster2.list <- filter(top500_combinedroc, cluster == 2)$gene
combined.cluster3.list <- filter(top500_combinedroc, cluster == 3)$gene
combined.cluster4.list <- filter(top500_combinedroc, cluster == 4)$gene
combined.cluster5.list <- filter(top500_combinedroc, cluster == 5)$gene

length(intersect(combined.cluster0.list, day1cluster0.list))
length(intersect(combined.cluster0.list, day1cluster2.list))
length(intersect(combined.cluster1.list, combined.cluster4.list))

day87cluster1.list <- filter(top500_87roc, cluster == 1)$gene
day87cluster2.list <- filter(top500_87roc, cluster == 2)$gene
day87cluster3.list <- filter(top500_87roc, cluster == 3)$gene
day87cluster5.list <- filter(top500_87roc, cluster == 5)$gene
length(intersect(day1cluster0.list, day1cluster2.list))

day87intersectionc2c4 <- intersect(day87cluster2.list, day87cluster4.list)
length(intersect(day87intersectionc2c4, day87cluster0.list))

venn.combined.df2 <- data.frame(cluster0 = c(500,0,12,6,15), cluster1 = c(0,500,31,87,54), 
  cluster3 = c(12,31,500,30,43), cluster4 = c(6,87,30,500,35), cluster5 = c(15,54,43,35,500))
row.names(venn.combined.df2) <- c(0,1,3,4,5)
View(venn.combined.df2)
heatmap(as.matrix(venn.combined.df2))
heatmap(as.matrix(venn.combined.df2), scale="column", col = terrain.colors(256))


day87cartilage.markers <- FindMarkers(sobj87, ident.1 = 3, ident.2 = c(0, 2, 4, 5), min.pct = 0.25)
day87cartilage.markers2 <- rownames_to_column(day87cartilage.markers, var = "geneid")
colnames(day87cartilage.markers2)[c(2, 3, 4, 5, 6)] <- paste('cartilage.day87', colnames(day87cartilage.markers2)[c(2, 3, 4, 5, 6)], sep = '_')
colnames(geneid.protein.nodedata.rowz)[c(2:26)] <- paste('TMT', colnames(geneid.protein.nodedata.rowz)[c(2:26)], sep = '_')
cartilage.day87.motherlode <- left_join(day87cartilage.markers2, geneid.protein.nodedata.rowz, by="geneid")
cartilage.day87.motherlodeREV <- left_join(geneid.protein.nodedata.rowz, day87cartilage.markers2, by="geneid")

# ggplot(filter(cartilage.day87.motherlode, cartilage.day87_avg_log2FC > 1)) +
ggplot(cartilage.day87.motherlode) +
          geom_point(aes(x = log10(TMT_NCSCmedians), y = log10(TMT_orgmedians), color = cartilage.day87_avg_log2FC), size = 5) +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          ylim(-4,4)

day87melanocytes.markers <- FindMarkers(sobj87, ident.1 = 5, ident.2 = c(0, 2, 3, 4), min.pct = 0.25)
day87melanocytes.markers2 <- rownames_to_column(day87melanocytes.markers, var = "geneid")
colnames(day87melanocytes.markers2)[c(2, 3, 4, 5, 6)] <- paste('melanocytes.day87', colnames(day87melanocytes.markers2)[c(2, 3, 4, 5, 6)], sep = '_')
melanocytes.day87.motherlode <- left_join(day87melanocytes.markers2, geneid.protein.nodedata.rowz, by="geneid")
melanocytes.day87.motherlodeREV <- left_join(geneid.protein.nodedata.rowz, day87melanocytes.markers2, by="geneid")

# ggplot(filter(melanocytes.day87.motherlode, melanocytes.day87_avg_log2FC > 1)) +
ggplot(melanocytes.day87.motherlode) +
          geom_point(aes(x = log10(TMT_NCSCmedians), y = log10(TMT_orgmedians), color = melanocytes.day87_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          ylim(-4,4)

day87neuroids.markers <- FindMarkers(sobj87, ident.1 = c(0, 2, 4), ident.2 = 3, min.pct = 0.25)
day87neuroids.markers2 <- rownames_to_column(day87neuroids.markers, var = "geneid")
colnames(day87neuroids.markers2)[c(2, 3, 4, 5, 6)] <- paste('neuroids.day87', colnames(day87neuroids.markers2)[c(2, 3, 4, 5, 6)], sep = '_')
neuroids.day87.motherlode <- left_join(day87neuroids.markers2, geneid.protein.nodedata.rowz, by="geneid")
neuroids.day87.motherlodeREV <- left_join(geneid.protein.nodedata.rowz, day87neuroids.markers2, by="geneid")

# ggplot(filter(neuroids.day87.motherlode, neuroids.day87_avg_log2FC > 1)) +
ggplot(neuroids.day87.motherlode) +
          geom_point(aes(x = log10(TMT_NCSCmedians), y = log10(TMT_orgmedians), color = neuroids.day87_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          ylim(-4,4)D

day1ncsc.markers <- FindMarkers(sobj.combined, ident.1 = 0, ident.2 = c(1,3,4,5), min.pct = 0.25)
day1ncsc.markers2 <- rownames_to_column(day1ncsc.markers, var = "geneid")
colnames(day1ncsc.markers2)[c(2, 3, 4, 5, 6)] <- paste('ncsc.day1', colnames(day1ncsc.markers2)[c(2, 3, 4, 5, 6)], sep = '_')
ncsc.day1.motherlode <- left_join(day1ncsc.markers2, geneid.protein.nodedata.rowz, by="geneid")
ncsc.day1.motherlodeREV <- left_join(geneid.protein.nodedata.rowz, day1ncsc.markers2, by="geneid")

################################### NEW PLOTTING COLLAGEN ##########################

ggplot(filter(ncsc.day1.motherlodeREV, grepl('^COL[0-9]{1,2}A[0-9]$', geneid) == TRUE)) +
          geom_point(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = ncsc.day1_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(axis.text.x = element_text(size = 24),
                axis.text.y = element_text(size = 23),
                axis.title.x = element_text(size = 24),
                axis.title.y = element_text(size = 24)) +
          ggtitle('ncsc.day1.motherlodeREV)') 

################################### NEW PLOTTING ###################################

# ggplot(filter(ncsc.day1.motherlodeREV, ncsc.day1_avg_log2FC > 1)) +
# ggplot(ncsc.day1.motherlodeREV) +
ggplot() +
          geom_point(data = filter(ncsc.day1.motherlodeREV, is.na(snRNAseq_avg_log2FC)==TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(ncsc.day1.motherlodeREV, is.na(snRNAseq_avg_log2FC)==FALSE & snRNAseq_avg_log2FC > 0.1),
          # geom_point(data = filter(ncsc.day1.motherlodeREV, is.na(snRNAseq_avg_log2FC)==FALSE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('ncsc.day1.motherlodeREV') +
          ylim(-1, 1) +
          xlim(-1, 1)

ggplot() +
          geom_point(data = filter(neuroids.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "gray", size = 2) +
          geom_point(data = filter(neuroids.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE & snRNAseq_avg_log2FC > 0.1),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('neuroids.day87.motherlodeREV') +
          ylim(-1, 1) +
          xlim(-1, 1)

ggplot() +
          geom_point(data = filter(melanocytes.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "gray", size = 2) +
          geom_point(data = filter(melanocytes.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE & snRNAseq_avg_log2FC > 0.1),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(data = filter(melanocytes.day87.motherlodeREV, (geneid %in% miniCOLmarkers) == TRUE), 
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('melanocytes.day87.motherlodeREV') +
          ylim(-1, 1) +
          xlim(-1, 1)

mini.elastic <-c("CDH11", "ELN", "COL2A1", "SOX9")
VlnPlot(sobj.combined, features = mini.elastic)
elastic.cartilage.markers <- c("ELN", "LUM", "BGN", "DCN", "VCAN", "FBN1", "LOXL1", "COL1A1", "COL2A1", "COL6A2", "COL6A3")
VlnPlot(sobj.combined, features = elastic.cartilage.markers)
# Warning message:
# In FetchData.Seurat(object = object, vars = features, slot = slot) :
#   The following requested variables were not found: FBLN3

ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE & snRNAseq_avg_log2FC > 0.1),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV')
          #  +
          # ylim(-1, 1) +
          # xlim(-1, 1)

ggplot() +
          geom_point(data = filter(ncsc.day1.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(ncsc.day1.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('ncsc.day1.motherlodeREV)') +
          ylim(-1, 1) +
          xlim(-1, 1)

#########################################################################
#### In vitro elastic cartilage reconstruction using human auricular perichondrial chondroprogenitor cell–derived micro 3D spheroids
#### Takayoshi et al. 2022. Journal of Tissue Engineering 13: 1-13.
#########################################################################

VlnPlot(sobj.combined, features = c("CD90", "CD44", "CD73", "CD105"))
# Warning message:
# In FetchData.Seurat(object = object, vars = features, slot = slot) :
#   The following requested variables were not found: CD90, CD73, CD105
VlnPlot(sobj1, features = c("THY1", "CD44", "NT5E", "ENG"))
# Warning message:
# In FetchData.Seurat(object = object, vars = features, slot = slot) :
#   The following requested variables were not found: CD90, CD73, CD105


##########################################################################
############# using Mark's allmarkers chondrocytes #######################
##########################################################################

allmarkers <- c("ACAN", "AGRN", "ASPN", "BGN", "COL11A1", "COL11A2", "COL12A1", "COL14A1", "COL15A1", "COL16A1", "COL18A1", "COL1A1", "COL1A2", "COL22A1", "COL24A1", "COL26A1", "COL2A1", "COL3A1", "COL4A1", "COL4A2", "COL4A3BP", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL9A1", "COLGALT1", "COLGALT2", "COMP", "CSPG4", "CTGF", "CTHRC1", "CYR61", "DAG1", "DCN", "ECM1", "EDIL3", "EFEMP1", "EFEMP2", "EMILIN1", "EMILIN3", "F9", "FAM20B", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FN1", "FRAS1", "FREM2", "HABP4", "HAPLN1", "HAPLN3", "HMCN1", "HMMR", "HSPG2", "LAMA1", "LAMA3", "LAMB1", "LAMB2", "LAMC1", "LGALS3BP", "LTBP1", "LTBP2", "LTBP3", "LTBP4", "LUM", "LYSMD1", "LYSMD2", "LYSMD3", "MATN1", "MATN2", "MATN3", "MATN4", "MCOLN1", "MCOLN3", "MFAP1", "MFAP2", "MFAP4", "MGP", "NID1", "NID2", "OBSL1", "OGN", "OLFML3", "PCOLCE", "PIGU", "PIGW", "PIGX", "PRG4", "SEMA3E", "SGCD", "SGCE", "SMOC1", "SPARC", "SPOCK1", "SPON2", "TFIP11", "THBS1", "THBS2", "THBS3", "THBS4", "TIMP1", "TIMP2", "TIMP3", "TNN", "TUFT1", "VCAN", "VTN", "VWA1", "ANXA6", "CD44", "CD151", "ITM2A", "FOXC1", "FOXC2", "SOX5", "SOX6", "SOX9", "CTSB", "CHADL", "CHAD", "CRTAC1", "DSPG3", "IBSP", "MIA", "OTOR", "SRPX", "ELN")
ggplot() +
          geom_point(data = filter(neuroids.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(neuroids.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          (geneid %in% allmarkers) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('neuroids.day87.motherlodeREV') +
          ylim(-30, 100)

##########################################################################
############# using auricular chondrocytes #######################
# Differences in Cartilage-Forming Capacity of Expanded Human Chondrocytes from Ear and Nose and Their Gene Expression Profiles
##########################################################################

auricular.chondrocytes <- c("PRG4", "BMP5", "FGL2", "MAB21L2", "DLX5", "MEIS1", "SEMA3E", "FAM107A", "CCKAR", "OGN", "MEIS2", "CYP24A1", "PI15", "HMCN1", "EVI1", "GREM2", "PTGS1", "PCSK1", "CHODL", "AGPAT9", "mohawk", "LRRC17", "NDP", "ELN", "FLRT2", "NPTX2", "CGNL1", "LRFN5", "NTNG1", "PODXL", "ELMO1", "FBN2", "GNG11", "KIAA0746", "SFRP1", "SAMD9", "COLEC12", "CRIP1", "NOVA1", "CLGN", "PTPRU", "THBD", "EIF1AX", "RPS6KA5", "ETV1", "ERG", "E26", "HSD11B1", "NGEF", "TIAM1", "MATN2", "NPTX1", "CACNB2", "SPRY2", "GALNT6", "DLX6", "G0S2", "PITPNC1")
ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          (geneid %in% auricular.chondrocytes) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV') +
          ylim(-30, 100)

##########################################################################
############# using growth plate chondrocytes ############################
# Genes uniquely expressed in human growth plate chondrocytes uncover a distinct regulatory network 2017
##########################################################################

growthplate.chondrocytes <- c("RELB", "SOX9", "NKX3-1", "RARG", "HMGA1", "MAFF", "FOSL1", "ATF3", "SNAI2", "PPARD", "CEBPD", "SMAD7", "LEF1", "EGR2", "NFATC1", "EGR4", "EGR3", "BCL6", "ERG", "ETS2")
ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          (geneid %in% growthplate.chondrocytes) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV') +
          ylim(-30, 100)

# > intersect(growthplate.chondrocytes, allmarkers)
# [1] "SOX9"

########################################################
################# elastifiber.markers ##################
########################################################

elastifiber.markers <- c("ELN", "FBN1", "FBN2", "MFAP4", "LTBP4", "LOX", "LOXL1", "LOXL2", "EMILIN1", "LTBP1", "LTBP2", "LTBP3", "DCN")
ggplot() +
          geom_point(data = filter(melanocytes.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          geom_point(data = filter(melanocytes.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          (geneid %in% elastifiber.markers) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          geom_label_repel(data = filter(melanocytes.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          (geneid %in% elastifiber.markers) == TRUE),
          aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('melanocytes.day87.motherlodeREV') +
          ylim(-5, 40) +
          xlim(-5, 5)

########################################################
####### redo FindMarkers with test.use = "roc" #########
########################################################

day1ncsc.markersroc <- FindMarkers(sobj.combined, ident.1 = 0, ident.2 = c(1,3,4,5), min.pct = 0.25, test.use = "roc")
day1ncsc.markersroc2 <- rownames_to_column(day1ncsc.markersroc, var = "geneid")
colnames(day1ncsc.markersroc2)[c(2, 3, 4, 5, 6)] <- paste('ncsc.day1', colnames(day1ncsc.markersroc2)[c(2, 3, 4, 5, 6)], sep = '_')
ncsc.day1.motherlodeROC <- left_join(geneid.protein.nodedata.rowz, day1ncsc.markersroc2, by = "geneid")
colnames(ncsc.day1.motherlodeROC)
colnames(ncsc.day1.motherlodeROC)[30] <- "snRNAseq_avg_log2FC"

ggplot() +
          geom_point(data = filter(ncsc.day1.motherlodeROC, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          # geom_point(data = filter(neuroids.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          # (geneid %in% allmarkers) == TRUE),
          geom_point(data = filter(ncsc.day1.motherlodeROC, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('ncsc.day1.motherlodeROC') +
          ylim(-30, 100)

day87neuroids.markersroc <- FindMarkers(sobj87, ident.1 = c(0, 2, 4), ident.2 = 3,  min.pct = 0.25, test.use = "roc")
day87neuroids.markersroc2 <- rownames_to_column(day87neuroids.markersroc, var = "geneid")
colnames(day87neuroids.markersroc2)[c(2, 3, 4, 5, 6)] <- paste('neuroids.day87', colnames(day87neuroids.markersroc2)[c(2, 3, 4, 5, 6)], sep = '_')
neuroids.day87.motherlodeROC <- left_join(geneid.protein.nodedata.rowz, day87neuroids.markersroc2, by = "geneid")
colnames(neuroids.day87.motherlodeROC)
colnames(neuroids.day87.motherlodeROC)[30] <- "snRNAseq_avg_log2FC"

ggplot() +
          geom_point(data = filter(neuroids.day87.motherlodeROC, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          # geom_point(data = filter(neuroids.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          # (geneid %in% allmarkers) == TRUE),
          geom_point(data = filter(neuroids.day87.motherlodeROC, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('neuroids.day87.motherlodeROC') +
          ylim(-30, 100)

day87cartilage.markersroc <- FindMarkers(sobj87, ident.1 = 3, ident.2 = c(0, 2, 4, 5),  min.pct = 0.25, test.use = "roc")
day87cartilage.markersroc2 <- rownames_to_column(day87cartilage.markersroc, var = "geneid")
colnames(day87cartilage.markersroc2)[c(2, 3, 4, 5, 6)] <- paste('cartilage.day87', colnames(day87cartilage.markersroc2)[c(2, 3, 4, 5, 6)], sep = '_')
cartilage.day87.motherlodeROC <- left_join(geneid.protein.nodedata.rowz, day87cartilage.markersroc2, by = "geneid")
colnames(cartilage.day87.motherlodeROC)
colnames(cartilage.day87.motherlodeROC)[30] <- "snRNAseq_avg_log2FC"

ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeROC, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "grey", size = 2) +
          # geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          # (geneid %in% allmarkers) == TRUE),
          geom_point(data = filter(cartilage.day87.motherlodeROC, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeROC') +
          ylim(-30, 100)

day87melanocytes.markersroc <- FindMarkers(sobj87, ident.1 = 5, ident.2 = c(0, 2, 3, 4),  min.pct = 0.25, test.use = "roc")
day87melanocytes.markersroc2 <- rownames_to_column(day87melanocytes.markersroc, var = "geneid")
colnames(day87melanocytes.markersroc2)[c(2, 3, 4, 5, 6)] <- paste('melanocytes.day87', colnames(day87melanocytes.markersroc2)[c(2, 3, 4, 5, 6)], sep = '_')
melanocytes.day87.motherlodeROC <- left_join(geneid.protein.nodedata.rowz, day87melanocytes.markersroc2, by = "geneid")
colnames(melanocytes.day87.motherlodeROC)
colnames(melanocytes.day87.motherlodeROC)[30] <- "snRNAseq_avg_log2FC"

ggplot() +
          geom_point(data = filter(melanocytes.day87.motherlodeROC, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians), color = "grey", size = 2) +
          # geom_point(data = filter(melanocytes.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
          # (geneid %in% allmarkers) == TRUE),
          geom_point(data = filter(melanocytes.day87.motherlodeROC, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = TMT_NCSCmedians, y = TMT_orgmedians, label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('melanocytes.day87.motherlodeROC') +
          ylim(-30, 100)

########################################################
########### ggplot with log2 transformation ############
########################################################

ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "grey", size = 2) +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE & snRNAseq_avg_log2FC > 0.1),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV') +
          xlim(-10, 5) +
          ylim(-10, 10)

########################################################
####### add contour map of geom_point ########
########################################################

ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "grey", size = 2) +
          geom_density_2d(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "red") +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), label = geneid), size = 3) +
          # scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV') +
          xlim(-10, 5) +
          ylim(-10, 10)

ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "grey", size = 2) +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), color = snRNAseq_avg_log2FC), size = 3) +
          # geom_density_2d(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
          #   aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "red") +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          # geom_label_repel(aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), label = geneid), size = 3) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          # theme(legend.position = 'None',
          theme(
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV') +
          xlim(-10, 5) +
          ylim(-10, 10)


########################################################
############### cell cycle stage scoring ###############
########################################################

# scoring the sequencing data relative to the cell cycle genes
CellCycleScoring(sobj.combined,
                  s.features = cc.genes.updated.2019$s.genes,
                  g2m.features = cc.genes.updated.2019$g2m.genes,
                  set.ident = TRUE) -> sobj.combined

# result
sobj.combined[[]]

# simple bar graph
as_tibble(sobj.combined[[]]) %>% ggplot(aes(Phase)) + geom_bar()
# scatter plot
as_tibble(sobj.combined[[]]) %>% ggplot(aes(x = S.Score, y = G2M.Score, color = Phase)) + geom_point()

# bar graph with clusters and raw numbers

as_tibble(sobj.combined[[]]) %>% ggplot(aes(x = seurat_clusters, fill = Phase)) +
                                       geom_bar()
# bar graph with clusters and reative proportion per cluster

as_tibble(sobj.combined[[]]) %>% ggplot(aes(x = seurat_clusters, fill = Phase)) +
                                       geom_bar(position = "fill", stat= "count")

# tables
table(sobj.combined@meta.data$Phase)
table(sobj.combined@meta.data$seurat_clusters)

# The following line resets the "default" view.
sobj.combined <- SetIdent(sobj.combined, value = sobj.combined@meta.data$seurat_clusters)

########################################################
####### fetch data ########
########################################################

# FetchData can pull anything from expression matrices, cell embeddings, or metadata
# FetchData(object = pbmc, vars = c("PC_1", "percent.mito", "MS4A1"))

FetchData(object = sobj.combined, vars = c("PC_1"))

PCAPlot(sobj.combined)
# is equivalent to
ggplot(FetchData(object = sobj.combined, vars = c("PC_1", "PC_2", "seurat_clusters"))) +
  geom_point(aes(x = PC_1, y = PC_2, color = seurat_clusters))
 
ggplot(FetchData(object = sobj.combined, vars = c("UMAP_1", "UMAP_2", "Phase"))) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = Phase))

ggplot(FetchData(object = sobj.combined, vars = c("UMAP_1", "UMAP_2", "seurat_clusters"))) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters))

# Embeddings(object = sobj.combined, reduction = "pca")
Loadings(object = sobj.combined, reduction = "pca")

FeatureScatter(object = sobj.combined, feature1 = "GREM1", feature2 = "PC_2")
RidgePlot(object = sobj.combined, feature = c("GREM1", "GREM2"))



specifmarkers1 <- setdiff(combined.cluster1.list, intersect(combined.cluster1.list,combined.cluster4.list))
specifmarkers1 <- setdiff(specifmarkers1, intersect(combined.cluster1.list,combined.cluster0.list))
specifmarkers1 <- setdiff(specifmarkers1, intersect(combined.cluster1.list,combined.cluster3.list))
specifmarkers1 <- setdiff(specifmarkers1, intersect(combined.cluster1.list,combined.cluster5.list))

specifmarkers3 <- setdiff(combined.cluster3.list, intersect(combined.cluster3.list,combined.cluster0.list))
specifmarkers3 <- setdiff(specifmarkers3, intersect(combined.cluster3.list,combined.cluster1.list))
specifmarkers3 <- setdiff(specifmarkers3, intersect(combined.cluster3.list,combined.cluster4.list))
specifmarkers3 <- setdiff(specifmarkers3, intersect(combined.cluster3.list,combined.cluster5.list))

specifmarkers5 <- setdiff(combined.cluster5.list, intersect(combined.cluster5.list,combined.cluster0.list))
specifmarkers5 <- setdiff(specifmarkers5, intersect(combined.cluster5.list,combined.cluster1.list))
specifmarkers5 <- setdiff(specifmarkers5, intersect(combined.cluster5.list,combined.cluster3.list))
specifmarkers5 <- setdiff(specifmarkers5, intersect(combined.cluster5.list,combined.cluster4.list))

########################################################
######### work on Friday April 14th, 2023 ####################
########################################################

filter(geneid.protein.nodedata.rowz, (geneid %in% c("GREM1", "GREM2", "DKK1", "KREMEN1", "KREMEN2")) == TRUE) %>% select(1, 22, 23)

library(ggrepel)
elastic.fiber.genes <- c("ELN", "FBN1", "FBN2", "MFAP4", "LOX", "LOXL1", "LOXL2",
                        "EMILIN1", "LTBP1", "LTBP2", "LTBP3", "LTBP4", "DCN")
ggplot() +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians)), color = "grey", size = 2) +
          geom_point(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
                                                                   (geneid %in% elastic.fiber.genes) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), color = snRNAseq_avg_log2FC), size = 4) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          geom_label_repel(data = filter(cartilage.day87.motherlodeREV, is.na(snRNAseq_avg_log2FC) == FALSE &
                                                                   (geneid %in% elastic.fiber.genes) == TRUE),
            aes(x = log2(TMT_NCSCmedians), y = log2(TMT_orgmedians), label = geneid), max.overlaps = Inf, size = 4) +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          theme_bw() +
          theme(axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('cartilage.day87.motherlodeREV') +
          xlim(-10, 5) +
          ylim(-10, 10)

########################################################
######### work on Monday April 17th, 2023 ##############
########################################################

library(devtools)
BiocManager::install("ComplexHeatmap")
install_github('msraredon/Connectome', ref = 'master')
library(Seurat)
library(tidyverse)
library(Connectome)
library(cowplot)

VlnPlot(sobj87, features = c("nCount_RNA", "nFeature_RNA"))
ggplot(FetchData(object = sobj87, vars = c("UMAP_1", "UMAP_2", "seurat_clusters"))) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters))

# panc8 <- NormalizeData(panc8)
# connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
# genes <- connectome.genes[connectome.genes %in% rownames(panc8)]
# panc8 <- ScaleData(panc8,features = genes)
# panc8.con <- CreateConnectome(panc8,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)

connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                        Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
# should I have redone the NormalizeData ? 
# I didn't.

congenes <- connectome.genes[connectome.genes %in% rownames(sobj87)]
sobj87 <- ScaleData(sobj87,features = congenes)
sobj87.con <- CreateConnectome(sobj87, species = 'human', min.cells.per.ident = 75,
                                      p.values = F, calculate.DOR = F)
p1 <- ggplot(sobj87.con, aes(x = ligand.scale)) + geom_density() + ggtitle('Ligand.scale')
p2 <- ggplot(sobj87.con, aes(x = recept.scale)) + geom_density() + ggtitle('Recept.scale')
p3 <- ggplot(sobj87.con, aes(x = percent.target)) + geom_density() + ggtitle('Percent.target')
p4 <- ggplot(sobj87.con, aes(x = percent.source)) + geom_density() + ggtitle('Percent.source')
plot_grid(p1, p2, p3, p4)

sobj87.con2 <- FilterConnectome(sobj87.con, min.pct = 0.1, min.z = 0.25, remove.na = T)
# Pre-filter edges:  60550
# Post-filter edges:  599
# Connectome filtration completed

install.packages("gridGraphics")
library(gridGraphics)
p1 <- NetworkPlot(sobj87.con2,features = 'PDGFC', min.pct = 0.1, weight.attribute = 'weight_sc', include.all.nodes = T)
p2 <- NetworkPlot(sobj87.con2,features = 'PDGFC', min.pct = 0.75,weight.attribute = 'weight_sc', include.all.nodes = T)
plot_grid(p1,p2,nrow=1)

p1 <- NetworkPlot(sobj872.con,features = 'TGFB', min.pct = 0.1, weight.attribute = 'weight_sc', include.all.nodes = T)
plot_grid(p1, nrow=1)

Centrality(sobj87.con2,
           modes.include = NULL,
           min.z = NULL,
           weight.attribute = 'weight_sc',
           group.by = 'mode')

Centrality(sobj87.con2,
           modes.include = c('PDGF'),
           weight.attribute = 'weight_sc',
           min.z = 0,
           group.by = 'mechanism')

# First, let’s select the top 10 signaling vectors for each cell-cell vector:

test <- sobj87.con2
test <- data.frame(test %>% group_by(vector) %>% top_n(10,weight_sc))

# Then, let’s focus in on specific cell types of interest:

cells.of.interest <- c('0', '1', '3', '4')

# Using edgeweights from normalized slot:
CircosPlot(test, weight.attribute = 'weight_norm',
            sources.include = cells.of.interest,
            targets.include = cells.of.interest,
            lab.cex = 0.6,
            title = 'Edgeweights from normalized slot')

# Using edgeweights from scaled slot:A
CircosPlot(test, weight.attribute = 'weight_sc',
            sources.include = cells.of.interest,
            targets.include = cells.of.interest,
            lab.cex = 0.6,
            title = 'Edgeweights from scaled slot')

# Using separate ligand and receptor expression (from normalized slot)
CircosPlot(test,weight.attribute = 'weight_norm',
            sources.include = cells.of.interest,
            targets.include = cells.of.interest,
            balanced.edges = F,
            lab.cex = 0.6,
            title = 'Ligand vs. receptor expression (from normalized slot)')
# PROBLEM?
# Pre-filter edges:  171
# Post-filter edges:  100
# Connectome filtration completed
# There are more than one numeric columns in the data frame.
# Take the first two numeric columns and draw the link ends
# with unequal width.
# Type `circos.par$message = FALSE` to suppress the message.
# Warning message:
# Since you have set `order`, you should better set `grid.col`
# as a named vector where sector names are the vector names,
# or else the color will be wrongly assigned. 

# Using separate ligand and receptor expression (from scaled slot)
CircosPlot(test,weight.attribute = 'weight_sc',
            sources.include = cells.of.interest,
            targets.include = cells.of.interest,
            balanced.edges = F,
            lab.cex = 0.6,
            title = 'Ligand vs. receptor expression (from scaled slot)')
# PROBLEM?
# Pre-filter edges:  171
# Post-filter edges:  100
# Connectome filtration completed
# There are more than one numeric columns in the data frame.
# Take the first two numeric columns and draw the link ends
# with unequal width.
# Type `circos.par$message = FALSE` to suppress the message.
# Warning message:
# Since you have set `order`, you should better set `grid.col`
# as a named vector where sector names are the vector names,
# or else the color will be wrongly assigned. 

CircosPlot(test, targets.include = '3', lab.cex = 0.6)

CircosPlot(test, sources.include = '3', lab.cex = 0.6)

########################################################
########################################################
########################################################
VizDimLoadings(sobj.combined, dims = 1:2, reduction = "pca")
# try plot_ly for 3D plots
plot_ly(x=temp, y=pressure, z=dtime, type="scatter3d", mode="markers", color=temp)

########################################################
####### make animated GIF with magick R package ########
########################################################

dataframe.list <- list(ncsc.day1.motherlodeREV, neuroids.day87.motherlodeREV,
              melanocytes.day87.motherlodeREV, cartilage.day87.motherlodeREV)

for (item in dataframe.list) {
p <- ggplot() +
          geom_point(data = filter(item, is.na(snRNAseq_avg_log2FC) == TRUE), 
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 2) +
          geom_point(data = filter(item, is.na(snRNAseq_avg_log2FC) == FALSE),
            aes(x = TMT_NCSCmedians, y = TMT_orgmedians, color = snRNAseq_avg_log2FC), size = 3) +
          geom_abline(slope = 1, linewidth = 1, color = "blue") +
          geom_hline(yintercept = 0, linewidth = 1, color = "orange") +
          geom_vline(xintercept = 0, linewidth = 1, color = "orange") +
          scale_color_gradientn(colours = rainbow(12), limits = c(-6, 6)) +
          ylim(-30, 100) +
          theme_bw() +
          theme(legend.position = 'None',
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18)) +
          ggtitle('XXXXXXXXXXX')
dfname <- paste0(item, ".png", sep = "")
ggsave(plot = p, filename = dfname, device = "png")
}


########################################
# min.pct
#     only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.
#     Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
# cluster6.markers <- FindMarkers(sobj87, ident.1 = 6, ident.2 = c(0, 1, 2, 3, 4, 5), min.pct = 0.25)
# cluster5.markers <- FindMarkers(sobj87, ident.1 = 5, ident.2 = c(0, 1, 2, 3, 4, 6), min.pct = 0.25)
# cluster4.markers <- FindMarkers(sobj87, ident.1 = 4, ident.2 = c(0, 1, 2, 3, 5, 6), min.pct = 0.25)
# cluster3.markers <- FindMarkers(sobj87, ident.1 = 3, ident.2 = c(0, 1, 2, 4, 5, 6), min.pct = 0.25)
# cluster2.markers <- FindMarkers(sobj87, ident.1 = 2, ident.2 = c(0, 1, 3, 4, 5, 6), min.pct = 0.25)
# cluster1.markers <- FindMarkers(sobj87, ident.1 = 1, ident.2 = c(0, 2, 3, 4, 5, 6), min.pct = 0.25)
# cluster0.markers <- FindMarkers(sobj87, ident.1 = 0, ident.2 = c(1, 2, 3, 4, 5, 6), min.pct = 0.25)



# FeaturePlot(sobj87, c("COL2A1", "COL9A1", "COL11A1", "COL11A2", "COL1A1", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "COL10A1"), max.cutoff = 'q95')
# VlnPlot(sobj87, c("COL2A1", "COL9A1", "COL11A1", "COL11A2", "COL1A1", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "COL10A1"))
# FeaturePlot(sobj87, c("PDGFRA", "PDGFRB", "PDGFA", "PDGFB", "PDGFC", "PDGFD"), max.cutoff = 'q95')
# VlnPlot(sobj87, c("PDGFRA", "PDGFRB", "PDGFA", "PDGFB", "PDGFC", "PDGFD"))

# library(SingleR)
# library(scater)
# library(loomR)
# library(celldex)
# library(scRNAseq)
# sobj87.sce <- as.SingleCellExperiment(sobj87)
# hpca.se <- HumanPrimaryCellAtlasData()
# hpca.se
# sobj87.hesc <- SingleR(test = sobj87.sce, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
# sobj87.hesc
# table(sobj87.hesc$labels)

#            Astrocyte Embryonic_stem_cells                  MSC Neuroepithelial_cell              Neurons          Osteoblasts 
#                   21                    3                    5                 3793                  131                    1 
#  Smooth_muscle_cells 
#                    1 

# VlnPlot(sobj87, c("COL1A1", "COL2A1", "ELN", "VCAN", "ACAN", "LUM", "BGN", "FBN1"))
# VlnPlot(sobj87, c("COL1A1", "COL2A1", "ELN", "LUM", "BGN", "FBN1"))

# VlnPlot(sobj87, features=c("PTN", "BDNF", "GDNF", "NGF", "VEGF", "MBP", "NLGN1", "PRKCI", "NEFL"))