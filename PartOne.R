#Seuralysis pilot  on env "dotsix" , with R3.6.1, Seurat v3, umap-learn, clustree
# Run it in /Seuralysis with :
#  echo "Rscript PartOne.R pilot ../outs" | qsub -V -cwd
# first argument is the run name for files, the second is the part of path between the main directory and the /filtered_feature_bc_matrix

args = commandArgs(trailingOnly=TRUE)
runam=args[1]
path=args[2]

runam
path

library(Seurat)
library(viridis)
library(ggplot2)
library(dplyr)
library(Matrix)
library(clustree)


pdf(file=paste("seur1_",runam,".pdf",sep=""))
pbmc_data <- Read10X(data.dir =paste("./../",path,"/filtered_feature_bc_matrix",sep=""))
pbmc <- CreateSeuratObject(counts = pbmc_data, min.cells = 3, min.features = 200, project = "pbmcs")
#14093 cells

#Filter cells with % mitochondrial genes > 10%
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.1)
pbmc <- subset(pbmc, subset = percent.mt < 10)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.1)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


# Load the HTO count matrix
matrix_dir = "./../outHTOer2/umi_count/"              
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
htoclass <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
colnames(htoclass) = barcode.names$V1
# If bacode names dont end with "-1" like the barcodes from RNA data, use " colnames(htoclass)<-unlist(head(lapply(barcode.names,paste,"-1", sep="")) )"
rownames(htoclass) = feature.names$V1
#only keep the common cell barcodes found from RNA and HTO data
joint_bcsb <- intersect(colnames(pbmc@assays$RNA@data),colnames(htoclass))
pbmc<-SubsetData(pbmc, cells=colnames(pbmc@assays$RNA@data[,joint_bcsb]))


htoclass <- as.matrix(htoclass[,joint_bcsb])
# Discard the unmapped hto barcodes
htoclass <- htoclass[setdiff(rownames(htoclass), c("unmapped")), ]
# renaming each HTO to avoid long names
rownames(htoclass)<-sub('[-][ATCG]+$','',rownames(htoclass))

pbmc[["HTO"]] <- CreateAssayObject(counts = htoclass)
# Normalize by Centred log ratio transformation
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc<-MULTIseqDemux(pbmc, assay = "HTO",autoThresh = T, maxiter = 5, qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)
pbmc[["HTO_classification"]]<-pbmc@meta.data$MULTI_classification
table(pbmc@meta.data$MULTI_ID)

classi<-c()
for(i in 1:length(pbmc@meta.data$MULTI_ID)){
	if (pbmc@meta.data$MULTI_ID[i]=="HTO1"|pbmc@meta.data$MULTI_ID[i]=="HTO4"|pbmc@meta.data$MULTI_ID[i]=="HTO7"|pbmc@meta.data$MULTI_ID[i]=="HTO8"|
	    pbmc@meta.data$MULTI_ID[i]=="HTO9"|pbmc@meta.data$MULTI_ID[i]=="HTO10"|pbmc@meta.data$MULTI_ID[i]=="HTO12"|pbmc@meta.data$MULTI_ID[i]=="HTO13"|
	    pbmc@meta.data$MULTI_ID[i]=="HTO14"|pbmc@meta.data$MULTI_ID[i]=="HTO15"){
		classi<-append(classi,"Singlet")
	} else if (pbmc@meta.data$MULTI_ID[i]=="Doublet"){
		classi<-append(classi,"Doublet")
	} else if (pbmc@meta.data$MULTI_ID[i]=="Negative"){
		classi<-append(classi,"Negative")
	}
}
pbmc[["HTO_classification.global"]]<-classi

# Plot the classification
plo<-HTOHeatmap(pbmc, assay = "HTO")
plo + scale_fill_viridis(name="Temperature")
#Keep only singlet cells
Idents(pbmc)<-pbmc@meta.data$HTO_classification.global
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
pbmc <- subset(pbmc, idents = "Singlet")

# pbmcd <- subset(pbmc, idents = "doublet")
# VlnPlot(pbmcd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.1)
# mean(pbmcd@meta.data$nCount_RNA)
# median(pbmcd@meta.data$nCount_RNA)
# mean(pbmc@meta.data$nCount_RNA)
# median(pbmc@meta.data$nCount_RNA)


#ADT assignation
matrix_dir= "./../outHA_Csc142/outADT/umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
pbmc_adt <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
colnames(pbmc_adt) = barcode.names$V1
rownames(pbmc_adt) = feature.names$V1
joint_bcsb <- intersect(colnames(pbmc@assays$RNA@data),colnames(pbmc_adt))
pbmc<-SubsetData(pbmc, cells=colnames(pbmc@assays$RNA@data[,joint_bcsb]))

pbmc_adt <- as.matrix(pbmc_adt[,joint_bcsb])
pbmc_adt <- pbmc_adt[setdiff(rownames(pbmc_adt), c("unmapped")), ]
rownames(pbmc_adt)<-sub('[-][ATCG]+$','',rownames(pbmc_adt))

pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc_adt)
pbmc<- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc<- ScaleData(pbmc, assay = "ADT")

# FeaturePlot(pbmc, features = rownames(pbmc_adt), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)


# Remove high read-count cells :
pbmc <- subset(pbmc, subset =nCount_RNA < (2*median(pbmc@meta.data$nCount_RNA)) & nFeature_RNA < (2*median(pbmc@meta.data$nFeature_RNA)))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc<- NormalizeData(pbmc)
# old way : (and was scaled only on variable features to be faster)
#pbmc_hashtagb <- FindVariableFeatures(pbmc, selection.method = "mean.var.plot")
#pbmc_hashtagb <- ScaleData(pbmc_hashtag, features = VariableFeatures(pbmc_hashtag))
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(pbmc), 10)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)
# Compute the number of PC to use (either explaining 95% variation or when 2 following PC have less than 0.1% difference of variation
elpt <-pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev)*100
cumul <- cumsum(elpt)
co1 <- which(cumul > 90 & elpt < 5)[1]
co2 <- sort(which((elpt[1:length(elpt) - 1] - elpt[2:length(elpt)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
if(pcs<10){pcs<-10}
pbmc <- FindNeighbors(pbmc, dims = 1:pcs)
pbmc<- RunUMAP(pbmc, dims = 1:pcs, umap.method= "umap-learn")
DimPlot(pbmc,reduction = "umap")
DefaultAssay(pbmc)<-"ADT"
FeaturePlot(pbmc,feature= rownames(pbmc_adt),sort=T)
DefaultAssay(pbmc)<-"RNA"
pbmcb <- FindClusters(pbmc, resolution =seq(from = 0.1, to =1.5, by = 0.1))


clustree(pbmcb,prefix="RNA_snn_res.")
clustree(pbmcb,prefix="RNA_snn_res.", node_colour="sc3_stability")

saveRDS(pbmc,file=paste("seur1_",runam,"_",pcs,"pcs.rds",sep=""))
dev.off()
