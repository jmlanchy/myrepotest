# from https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html
# Paragraph: Cell Cycle Scoring

# NOTE 1: I assume the data from mg-snrna-march2023 sequencing is already in the Rproject's Environment as a S4 object of class Seurat and named sobj.
# NOTE 2: "Seurat comes with a bunch of marker genes for different cell cycle stages"

library(Seurat)
library(tidyverse)
# test modif
# A look at the cell cycle marker genes
cc.genes.updated.2019

# scoring the sequencing data relative to the cell cycle genes
CellCycleScoring(sobj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> sobj

# result
sobj[[]]

# simple bar graph
as_tibble(sobj[[]]) %>% ggplot(aes(Phase)) + geom_bar()
# scatter plot
as_tibble(sobj[[]]) %>% ggplot(aes(x = S.Score, y = G2M.Score, color = Phase)) + geom_point()

# bar graph with clusters and raw numbers
as_tibble(sobj[[]]) %>% ggplot(aes(x = seurat_clusters, fill = Phase)) + geom_bar()
# bar graph with clusters and reative proportion per cluster
as_tibble(sobj[[]]) %>% ggplot(aes(x = seurat_clusters, fill = Phase)) + geom_bar(position = "fill", stat= "count")

# tables
table(sobj@meta.data$Phase)
table(sobj@meta.data$seurat_clusters)

# The following line resets the "default" view on the clusters (like UMAPPlot, TSNEPlot).
sobj <- SetIdent(sobj, value = sobj@meta.data$seurat_clusters)