library(stringr)
wd = getwd()
sample <- basename(normalizePath(wd))
phase = 'EFM'
    
library(getopt)
arg <- matrix(c("rawMatrix","r","1","character","raw matrix",
                "filterMatrix","f","1","character","filter matrix",
                "minCG","m","1","integer","minimum gene number for each cell, i.e. nFeature_RNA, default=0",
                "tfidfMin","t","1","integer","default = 1, Minimum value of tfidf to accept for a marker gene",
                "outdir","o","1","character","outdir of matrix"),byrow=T,ncol=5) 
opt = getopt(arg) 
if (is.null(opt$minCG)){
  opt$minCG <- 0
}            
if (is.null(opt$tfidfMin)){
  opt$tfidfMin <- 1
}
if (is.null(opt$rawMatrix)){
  opt$rawMatrix <- paste0("/data/input/",sample,"/Copy-scRNA-seq_v3.1.5/02.cDNAAnno/RawMatrix")
} 
if (is.null(opt$filterMatrix)){
  opt$filterMatrix <- paste0("/data/input/",sample,"/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix")
} 
if (is.null(opt$outdir)){
  opt$outdir <- paste0("/data/work/output/",sample,"/soupX_matrix/")
} 

library(DropletUtils)
library(SoupX)
library(Seurat)

options(future.globals.maxSize = 100000 * 1024^3)
# setwd("outdir1") #这样写会显示can not change work directory

toc <- Read10X(opt$filterMatrix,gene.column=1) 
tod <- Read10X(opt$rawMatrix,gene.column=1) 

#tod <- tod[rownames(toc),]

all <- toc
all <- CreateSeuratObject(all)
all <- subset(all, subset = nFeature_RNA > opt$minCG )
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:30)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:30)

matx <- all@meta.data
toc <- all@assays[["RNA"]]@counts

raw <- CreateSeuratObject(tod)
tod <- raw@assays[["RNA"]]@counts
tod <- tod[rownames(all),]

sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc, tfidfMin = opt$tfidfMin, forceAccept = T)

if(unique(sc$metaData$rho) > 0.5){
    sc = setContaminationFraction(sc, opt$rho)
}
samp <- c("Bone_marrow","PBMC","Pineal_gland")
if (length(grep("Bone_marrow|PBMC|Pineal_gland","$name"))==0) {
    if (unique(sc$metaData$rho) < 0.2) {
       sc = setContaminationFraction(sc, 0.2)
    }
}

out = adjustCounts(sc)

saveRDS(sc,paste0("/data/work/output/",sample,"/",sample,".rds"))

DropletUtils:::write10xCounts(opt$outdir, out,version="3")
