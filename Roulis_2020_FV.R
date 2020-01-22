## Drop-seq analysis

# Author: Rihao Qu/Jun Zhao
# Richard Flavell/Yuval Kluger Lab


library(Seurat)
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R")
source("./combineData.R")

# samples
sample_names <- c("F1A", "F1B", "F2", "V1A", "V1B", "V2")
n_sample <- length(sample_names)

# read in data
data_list <- vector("list", n_sample)
names(data_list) <- sample_names


for(i in 1:n_sample){
  data_list[[i]] <- read.table(
    paste(sample_names[i], "/gene_exon_tagged_dge.txt.gz", sep = ""), 
    sep = "\t", row.names = 1, header = T, stringsAsFactors = F
  )
}

data_exp <- combineData(
  data1 = data_list[[1]], data2 = data_list[[2]], name1 = sample_names[1], name2 = sample_names[2]
)

for(i in 3:n_sample){
  data_exp <- combineData(
    data1 = data_exp, data2 = data_list[[i]], name2 = sample_names[i]
  )
}

data_S <- CreateSeuratObject(raw.data = data_exp, min.genes = 0, names.delim = "_", names.field = 2)
strsplit2 <- function(x){strsplit(x,"_")[[1]][2]}
strsplit3 <- function(x){strsplit(x,"_")[[1]][1]}
data_S@meta.data$orig.ident <- sapply(rownames(data_S@meta.data), strsplit2)


# filter cells on nUMI
data_S <- FilterCells(data_S, subset.names = "nUMI", low.thresholds = 999)
table(data_S@meta.data$orig.ident)

# filter cells on mito.ratio
mito.genes <- grep(pattern = "mt-", x = rownames(x = data_S@data), value = TRUE, ignore.case = T)
data_S@meta.data$mito.ratio <- colSums(data_S@data[mito.genes, ])/colSums(data_S@data)
max(data_S@meta.data$mito.ratio)
VlnPlot(data_S, "mito.ratio", group.by = "orig.ident", point.size.use = 0.1)

data_S <- FilterCells(data_S, subset.names = "mito.ratio", high.thresholds = 0.1)
table(data_S@meta.data$orig.ident)

data_S@meta.data$cell_barcode <- sapply(rownames(data_S@meta.data), strsplit3)
barcode_F <- unique(data_S@meta.data$cell_barcode[is.element(data_S@meta.data$orig.ident,c("F1A","F1B","F2"))])
barcode_V <- unique(data_S@meta.data$cell_barcode[is.element(data_S@meta.data$orig.ident,c("V1A","V1B","V2"))])
write.table(barcode_F,"barcode_F.txt",quote=F,col.names=F,row.names=F)
write.table(barcode_V,"barcode_V.txt",quote=F,col.names=F,row.names=F)



# basic analysis
data_S <- NormalizeData(data_S)
data_S <- alraSeurat2(data_S)
data_S <- ScaleData(data_S, do.scale = F, do.center = T)
data_S <- RunPCA(data_S, pcs.compute = 20, pc.genes=rownames(data_S@data),do.print = F)
data_S <- RunTSNE(
  data_S, dims.use = 1:20, perplexity=50, check_duplicates = FALSE, do.fast = TRUE, seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/home/jz437/Tools/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
write.table(sub_data_S@dr$tsne@cell.embeddings[which(sub_data_S@meta.data$combined.ident=="F"),], "F_tsne.txt")
write.table(sub_data_S@dr$tsne@cell.embeddings[which(sub_data_S@meta.data$combined.ident=="V"),], "V_tsne.txt")


TSNEPlot(data_S, group.by = "orig.ident", pt.size = 1)
relabel <- function(x){
  if (is.element(x,c("F1A","F1B","F2"))) "F"
  else "V"
}
data_S@meta.data$combined.ident <- sapply(data_S@meta.data$orig.ident, relabel)
pdf("tsne.pdf")
TSNEPlot(data_S, group.by = "combined.ident", pt.size = 1.25, colors.use = c("blue", "red"))
dev.off()


# clustering
data_S <- FindClusters(data_S, reduction.type = "pca", dims.use = 1:36, resolution = 1, print.output = F, force.recalc = T)

pdf("cluster.pdf")
TSNEPlot(data_S, do.label = T, label.size = 5, pt.size = 1, colors.use = colors.use)
dev.off()
write.table(sub_data_S@meta.data$res.0.5[which(sub_data_S@meta.data$combined.ident=="F")], "F_clust_label.txt")
write.table(sub_data_S@meta.data$res.0.5[which(sub_data_S@meta.data$combined.ident=="V")], "V_clust_label.txt")

# cluster markers
nClust <- length(unique(data_S@ident))
cluster.markers <- FindAllMarkers(data_S, min.pct = 0.1, min.diff.pct = 0.09, only.pos = T)

