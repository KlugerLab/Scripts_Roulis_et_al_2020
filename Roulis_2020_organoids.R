## Drop-seq analysis

# Author: Jun Zhao
# Richard Flavell/Yuval Kluger Lab


library(Seurat)
library(dbscan)
library(RANN)
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R")
source("./combineData.R")


## Data preprocessing

# read data
drop_exp <- combineData(data1 = read.table("./DS014B_C/gene_exon_tagged_dge.txt.gz", 
                                           header = T,sep = "\t",stringsAsFactors = F, row.names = 1), 
                        data2 = read.table("./DS014B_T/gene_exon_tagged_dge.txt.gz",
                                           header = T,sep = "\t",stringsAsFactors = F, row.names = 1),
                        name1 = "C", name2 = "T")


# S obj
data_S <- CreateSeuratObject(raw.data = drop_exp, project = "organoid",
                             min.cells = 0, min.genes = 0, names.field = 2, names.delim = "_")
table(data_S@meta.data$orig.ident)


# filter cells
sum(data_S@meta.data$nUMI >= 1000)
data_S <- FilterCells(data_S, "nUMI", low.thresholds = 999)

mito.genes <- grep(pattern = "^mt-", x = rownames(x = data_S@data), value = TRUE, ignore.case = T)
data_S@meta.data$mito.ratio <- colSums(data_S@raw.data[mito.genes,data_S@cell.names])/
  colSums(data_S@raw.data[,data_S@cell.names])
max(data_S@meta.data$mito.ratio)
VlnPlot(data_S, "mito.ratio", group.by = "orig.ident", point.size.use = 0.1)

data_S <- FilterCells(data_S, subset.names = "mito.ratio", high.thresholds = 0.1)
VlnPlot(data_S, "mito.ratio", group.by = "orig.ident", point.size.use = 0.1)
table(data_S@meta.data$orig.ident)
length(data_S@cell.names)


# norm
data_S <- NormalizeData(data_S)

# ALRA
data_S <- alraSeurat2(data_S)

# analysis
data_S <- ScaleData(data_S, do.scale = F, do.center = T)
data_S <- RunPCA(data_S, pcs.compute = 29, pc.genes = rownames(data_S@data), genes.print = 5)

data_S <- RunTSNE(data_S, dims.use = 1:29, seed.use = 1)
TSNEPlot(data_S, group.by = "orig.ident", pt.size = 0.5)

data_S <- FindClusters(data_S, reduction.type = "pca", dims.use = 1:29, resolution = 0.8, 
                       print.output = F, force.recalc = T)
TSNEPlot(data_S, do.label = T, pt.size = 0.5)

# identify epithelial clusters
FeaturePlot(data_S, "Epcam", cols.use = c("blue","red"))



## Analyze epithelial cells

# get epi cells
ds014.epi_S <- FilterCells(data_S, subset.names = "nUMI",
                           cells.use = data_S@cell.names[!data_S@ident %in% c(0,6,9)])
length(ds014.epi_S@cell.names)

# redo analysis on epi cells
ds014.epi_S <- ScaleData(ds014.epi_S, do.scale = F, do.center = T)
ds014.epi_S <- RunPCA(ds014.epi_S, pcs.compute = 29, pc.genes = rownames(ds014.epi_S@data), genes.print = 5)

ds014.epi_S <- RunTSNE(ds014.epi_S, dims.use = 1:29, seed.use = 1)
TSNEPlot(ds014.epi_S, group.by = "orig.ident", pt.size = 0.5,
         colors.use = c("C" = "#CC0000", "T" = "#0000CC")) + 
  theme_cowplot() + theme(legend.title = element_blank())

ds014.epi_S <- FindClusters(ds014.epi_S, reduction.type = "pca", dims.use = 1:29, resolution = 1.8, 
                            print.output = F, force.recalc = T)
TSNEPlot(ds014.epi_S, do.label = T, pt.size = 0.5)



## DBSCAN on epi cells

# dbscan package
ds014.epi_dbscanObj <- dbscan(
  x = ds014.epi_S@dr$tsne@cell.embeddings, 
  eps = 2, minPts = 10
)
table(ds014.epi_dbscanObj$cluster)

ds014.epi_S@meta.data$dbscan <- ds014.epi_dbscanObj$cluster
TSNEPlot(ds014.epi_S, group.by = "dbscan", do.label = T)


# assign noisy points
for(i in 1:length(ds014.epi_S@cell.names)){
  if(ds014.epi_dbscanObj$cluster[i] == 0){
    nnObj <- nn2(data = ds014.epi_S@dr$tsne@cell.embeddings,
                 query = t(ds014.epi_S@dr$tsne@cell.embeddings[i,]),
                 k = 15)
    nnCluster <- ds014.epi_dbscanObj$cluster[nnObj$nn.idx[1,1:15]]
    ds014.epi_S@meta.data$dbscan[i] <- nnCluster[nnCluster != 0][1]
  }
}
table(ds014.epi_S@meta.data$dbscan)
TSNEPlot(ds014.epi_S, group.by = "dbscan", do.label = T)
nClust_epi <- length(unique(ds014.epi_S@meta.data$dbscan))

# merge clusters
ds014.epi_S@meta.data$dbscan.merge <- ds014.epi_S@meta.data$dbscan
ds014.epi_dbscan.order <- list(c(1,2,13), c(5,14,15,16), c(4,8), 9, 7, 6, 10, 11, 3, 12)

nClust_ds014.epi.merge <- length(ds014.epi_dbscan.order)
for(i in 1:nClust_ds014.epi.merge){
  ds014.epi_S@meta.data$dbscan.merge[ds014.epi_S@meta.data$dbscan %in% ds014.epi_dbscan.order[[i]]] <- i-1
}
table(ds014.epi_S@meta.data$dbscan.merge)
TSNEPlot(
  ds014.epi_S, group.by = "dbscan.merge", do.label = T, pt.size = 0.5, return.plot = T,
  colors.use = c("0"="#CC9900", "1"="#F8766D", "2"="#003399", "3"="#99CCFF", "4"="#669933",
                 "5"="#CCFF99", "6"="#9966CC", "7"="#9590FF", "8"="#990000", "9"="#FF99CC")
) + 
  theme_cowplot() + theme(legend.title = element_blank())


