## Drop-seq analysis

# Author: Jun Zhao
# Richard Flavell/Yuval Kluger Lab


library(Seurat)

source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R")
source("./combineData.R")


## Prepare data

# samples
sample_names <- c("C1A","C1B","C2A","C2B","C2C")
n_sample <- length(sample_names)

# read in data
data_list <- vector("list", n_sample)
names(data_list) <- sample_names

for(i in 1:n_sample){
  data_list[[i]] <- read.table(
    paste("./data/", sample_names[i], "/gene_exon_tagged_dge.txt.gz", sep = ""), 
    sep = "\t", row.names = 1, header = T, stringsAsFactors = F
  )
}



## Combine data directly

data_exp <- combineData(
  data1 = data_list[[1]], data2 = data_list[[2]], name1 = sample_names[1], name2 = sample_names[2]
)

for(i in 3:n_sample){
  data_exp <- combineData(
    data1 = data_exp, data2 = data_list[[i]], name2 = sample_names[i]
  )
}


# Seurat
data_S <- CreateSeuratObject(raw.data = data_exp, min.genes = 0, names.delim = "_", names.field = 2)

# transcript sum filtering
data_S <- FilterCells(data_S, subset.names = "nUMI", low.thresholds = 999)

# filter cells on mito.ratio
mito.genes <- grep(pattern = "mt-", x = rownames(x = data_S@data), value = TRUE, ignore.case = T)
data_S@meta.data$mito.ratio <- colSums(data_S@raw.data[mito.genes,data_S@cell.names])/
  colSums(data_S@raw.data[,data_S@cell.names])
max(data_S@meta.data$mito.ratio)
VlnPlot(data_S, "mito.ratio", group.by = "orig.ident", point.size.use = 0.1)

data_S <- FilterCells(data_S, subset.names = "mito.ratio", high.thresholds = 0.1)



## Analysis

# basic analysis
data_S <- NormalizeData(data_S)
data_S <- alraSeurat2(data_S)
data_S <- ScaleData(data_S, do.scale = F, do.center = T)

data_S <- RunPCA(data_S, pcs.compute = 36, pc.genes = rownames(data_S@data), do.print = F)
data_S <- RunTSNE(
  data_S, dims.use = 1:36, seed.use = 1
)
TSNEPlot(data_S, group.by = "orig.ident", pt.size = 0.5)

# clustering
data_S <- FindClusters(data_S, reduction.type = "pca", dims.use = 1:36, resolution = 1, print.output = F, force.recalc = T)
TSNEPlot(data_S, do.label = T, pt.size = 0.5)



## Remove immune cells

# check Ptprc
FeaturePlot(data_S, features.plot = "Ptprc", cols.use = rainbow(12, s = 0.7, alpha = .8)[9:1], pt.size = 0.5)

# remove immune clusters
data_mes_S <- SubsetData(data_S, cells.use = data_S@cell.names[!data_S@ident %in% c(4,6,10,12,14,16)], subset.raw = T)

data_mes_S <- ScaleData(data_mes_S, do.scale = F, do.center = T)
data_mes_S <- RunPCA(data_mes_S, pcs.compute = 36, pc.genes = rownames(data_mes_S@data), do.print = F)
plot(data_mes_S@dr$pca@sdev)

data_mes_S <- RunTSNE(
  data_mes_S, dims.use = 1:28, seed.use = 1
)

# clustering
data_mes_S <- FindClusters(data_mes_S, reduction.type = "pca", dims.use = 1:28, 
                           print.output = F, resolution = 1, force.recalc = T)
TSNEPlot(data_mes_S, do.label = T, pt.size = 0.5)
nClust_mes <- length(unique(data_mes_S@ident))


# reorder clusters
data_mes_S@meta.data$cluster <- c(
  "5"="0", "4"="1", "9"="2", "12"="3", "14"="4", "13"="5", "11"="6", "10"="7", "3"="8", 
  "2"="9", "1"="10", "0"="11", "7"="12", "6"="13", "8"="14"
)[data_mes_S@meta.data$res.1]
data_mes_S@meta.data$cluster <- factor(data_mes_S@meta.data$cluster, levels = as.character(c(0:14)))

cols_use <- c("#FD61D1","#E76BF3","#00BA38","#E58700","#A3A500","#00C0AF","#00BCD8","#619CFF",
              "#B983FF","#C99800","#00BF7D","#F8766D","#6BB100","#00B0F6","#FF67A4")
TSNEPlot(data_mes_S, do.label = T, pt.size = 0.5, group.by = "cluster", colors.use = cols_use) + 
  theme_cowplot() + theme(legend.title = element_blank())



## Ptgs2 clusters

# get Ptgs2+ cells
data_Ptgs2_S <- SubsetData(
  data_mes_S, cells.use = data_mes_S@cell.names[which(data_mes_S@data["Ptgs2",] > 0)], subset.raw = T
)


# analysis
data_Ptgs2_S <- ScaleData(data_Ptgs2_S, do.scale = F, do.center = T)
data_Ptgs2_S <- RunPCA(data_Ptgs2_S, pcs.compute = 36, pc.genes = rownames(data_Ptgs2_S@data), do.print = F)
plot(data_Ptgs2_S@dr$pca@sdev)

data_Ptgs2_S <- RunTSNE(
  data_Ptgs2_S, dims.use = 1:27, seed.use = 1
)
TSNEPlot(data_Ptgs2_S, group.by = "orig.ident", pt.size = 0.5)

# clustering
data_Ptgs2_S <- FindClusters(data_Ptgs2_S, reduction.type = "pca", dims.use = 1:27, 
                             print.output = F, force.recalc = T)
TSNEPlot(data_Ptgs2_S, do.label = T, pt.size = 0.5)
nClust_Ptgs2 <- length(unique(data_Ptgs2_S@ident))


# reorder cluster
data_Ptgs2_S@meta.data$cluster <- c(
  "5"="0", "6"="1", "1"="2", "3"="3", "0"="4", "2"="5", "4"="6"
)[data_Ptgs2_S@meta.data$res.0.8]
TSNEPlot(data_Ptgs2_S, do.label = T, pt.size = 0.5, group.by = "cluster") + 
  theme_cowplot() + theme(legend.title = element_blank())



