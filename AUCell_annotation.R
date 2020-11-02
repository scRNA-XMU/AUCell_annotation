#-----------------------------------------parameter-----------------------------------------
#input_path:存放cellranger 的3个输出文件（barcode，gene，matrix），包括一个sample信息的文件
#seurat_obj：生成的一个seurat对象
#nFeature_RNA：每个每个细胞表达的基因数
#nCount_RNA：每个细胞的UMI
#percent.mt：线粒体表达量
#normalization.method = "LogNormalize",scale.factor = 10000 ：  默认参数
#selection.method：'vst';'mean.var.plot';'dispersion '
#nfeatures:要筛选的高表达基因个数
#dispersion.cutoff：用于区分高表达的离散系数
#dims:是根据  ElbowPlot(seurat_obj)  这一步的结果设定的参数
#resolution：设置聚类的分辨率
#CANCER,B_CELL.....   :都是人为定义的marker基因，需要客服输入
#nCores：设定跑程序的CPU
#AUCannotations每个细胞的注释信息
#label.size：标注信息的字体大小
#highlights:需要高亮显示的细胞类型
#-----------------------------------------parameter-----------------------------------------





#-----------------------------------------设置工作路径-----------------------------------------
###设置R语言工作环境，用于存放输入数据和输出数据
dir = "/cluster/huanglab/cliang/lung_E_MTAB/"
setwd(dir)
#-----------------------------------------设置工作路径-----------------------------------------




#----------------------------------------相关包的载入----------------------------------------
library(AUCell)
library(Seurat)
if (!requireNamespace("BiocManager", quietly=TRUE))                            
  install.packages("BiocManager")                                              
# To support paralell execution:                                               
BiocManager::install(c("doMC", "doRNG"))                                       
n# For the main example:                                                       
BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"))        
# For the examples in the follow-up section of the tutorial:                   
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh",  
                       "dynamicTreeCut","R2HTML","Rtsne", "zoo"))              
#----------------------------------------相关包的载入----------------------------------------










#-----------------------------------载入相关数据，并进行seurat的标准处理流程----------------------------------------
input_path = "/cluster/huanglab/cliang/lung_E_MTAB/"
##这里不知道你们那边fastq生成count的时候时候可以生成一个这样的含有样本名的'sample_info.txt'， 将'sample_info.txt'读进来命名sample_info的矩阵
sample_info <- read.table(paste0(input_path,'sample_info.txt'), header = T, sep = "\t", stringsAsFactors = F)
##创一个list，将所有创建的seurat对象存在这里
objlist <- list()
for (i in 1:length(sample_info$SampleName)) {
  objlist[[i]] <- CreateSeuratObject(counts = Read10X(paste0(input_path, sample_info$SampleName[i], "outs/filtered_feature_bc_matrix")), 
                                     project = sample_info$SourceName[i])
}
##?????кܶ???count matrix??????Ҫ?????ľ??????ں???һ??
seurat_obj <- merge(x = objlist[[1]],
                    y = c(objlist[[2]],objlist[[3]]),
                    project = "lung")
##唯一的代码，只需输入seurat_obj
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
##理想是输出三张violin图，供客户选择阈值，这些nCount_RNA；nFeature_RNA；percent.mt的阈值需要输入
seurat_obj <- subset(seurat_obj,subset = nCount_RNA > 200 & nFeature_RNA < 6000 & nFeature_RNA > 100 & percent.mt < 10)  #理想输出一张图，供客户选择阈值
##normalization.method；scale.factor是需要输入的值，
seurat_obj <- NormalizeData(seurat_obj,normalization.method = "LogNormalize",scale.factor = 10000)  
###selection.method；mean.cutoff；dispersion.cutoff是需要输入的值
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mvp", 
                                   mean.cutoff = c(0.125, 3), dispersion.cutoff = c(0.5,1))
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj,features = VariableFeatures(object = seurat_obj))
####交互，输出一张图给客户，客户决定dims
ElbowPlot(seurat_obj) 
###客户输入dims
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:8)
###设置分辨率  resolution  需要输入
seurat_obj <- FindClusters(seurat_obj, resolution = 0.15)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:8)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:8)


cluster_graph_tsne <- DimPlot(seurat_obj, reduction = "tsne")
cluster_graph_umap <- DimPlot(seurat_obj, reduction = 'umap')
#-----------------------------------载入相关数据，并进行seurat的标准处理流程----------------------------------------





#-----------------------------------   AUCell包注释过程   ----------------------------------
#-----------------------------------0.加载scRNA-seq数据集和基因集----------------------------------
#AUCell的输入数据是表达矩阵和基因集
##这一步需要客户提供，或者我们有我们数据库

ALVEOLAR = c("CLDN18", "FOLR1", "AQP4", "PEBP4")
BCELL = c("CD79A", "IGKC", "IGLC3", "IGHG3")
EC = c("CLDN5", "FLT1", "CDH5", "RAMP2")
EPITHELIAL = c("CAPS", "TMEM190", "PIFO", "SNTN")
FIBROBLAST = c("COL1A1", "DCN", "COL1A2", "C1R")
MYELOID = c("LYZ", "MARCO", "CD68", "FCGR3A")
TCELL = c("CD3D", "TRBC1", "TRBC2", "TRAC")
CANCER = c("EPCAM", "ALCAM", "CD44", "PROM1")


ALVEOLAR.geneset = GeneSet(ALVEOLAR, setName="ALVEOLAR")
BCELL.geneset = GeneSet(BCELL, setName="BCELL")
EC.geneset = GeneSet(EC, setName="EC")
EPITHELIAL.geneset = GeneSet(EPITHELIAL, setName="EPITHELIAL")
FIBROBLAST.geneset = GeneSet(FIBROBLAST, setName="FIBROBLAST")
MYELOID.geneset = GeneSet(MYELOID, setName="MYELOID")
TCELL.geneset = GeneSet(TCELL, setName="TCELL")
CANCER.geneset = GeneSet(CANCER, setName="CANCER")
geneSets = GeneSetCollection(c(ALVEOLAR.geneset, BCELL.geneset, EC.geneset, EPITHELIAL.geneset, 
                               FIBROBLAST.geneset, MYELOID.geneset, TCELL.geneset, CANCER.geneset))
#-----------------------------------0.加载scRNA-seq数据集和基因集----------------------------------



#-----------------------------------   AUCell包注释过程   ----------------------------------
#----------------------------------确定具有给定基因特征或活性基因集的细胞----------------------------------
cells_rankings <- AUCell_buildRankings(seurat_obj[['RNA']]@counts, nCores=1, plotStats=TRUE)        
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)                                         
dataauc <- getAUC(cells_AUC)
#Determine the cells by max AUC
AUCannotations = apply(dataauc, 2, function(x){names(which.max(x))})
seurat_obj[["AUCannotations"]] = AUCannotations
Idents(seurat_obj) <- seurat_obj@meta.data$AUCannotations
####reduction  客户需要输入，想绘制的图，tsne或者umap   label：标签   label.size：字体大小
DimPlot(seurat_obj, reduction = "umap", label = T, label.size = 4)
###客户设置需要查看的细胞类型
highlights <- c("EC")
DimPlot(seurat_obj, reduction = "tsne", label = T, cells.highlight = colnames(seurat_obj)[seurat_obj@meta.data$AUCannotations==highlights])
#----------------------------------确定具有给定基因特征或活性基因集的细胞----------------------------------
#-----------------------------------   AUCell包注释过程   ----------------------------------


new.cluster.ids <- c('ALVEOLAR',
                     'BCELL',
                     'EC' ,
                     'EPITHELIAL' ,
                     'FIBROBLAST' ,
                     'MYELOID',
                     'TCELL',
                     'CANCER',)
names(new.cluster.ids) <- levels(seurat_obj@meta.data$seurat_clusters)
Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = 1) 
