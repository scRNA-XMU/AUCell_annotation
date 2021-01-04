#CANCER,B_CELL.....   :都是人为定义的marker基因，需要客服输入
#nCores：设定跑程序的CPU
#AUCannotations每个细胞的注释信息
#label.size：标注信息的字体大小
#highlights:需要高亮显示的细胞类型


#----------------------------------------相关包的载入----------------------------------------
library(AUCell)
library(Seurat)
#----------------------------------------相关包的载入----------------------------------------
#-----------------------------------------设置工作路径-----------------------------------------
###设置R语言工作环境，用于存放输入数据和输出数据
dir = "/cluster/huanglab/cliang/lung_E_MTAB/"
setwd(dir)
#-----------------------------------------设置工作路径-----------------------------------------




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
