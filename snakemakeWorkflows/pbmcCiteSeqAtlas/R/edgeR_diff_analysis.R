
library(edgeR)
library(getopt)
library(dplyr)
library(Seurat)
print(getwd())
source("figures/fun_r/seurat2PB.R")

spec <- matrix(c(
  'citeSeq_data_path',  'c', 1, "character",
  'output_path',  'o', 1, "character",
  'help',   'h', 0, "logical"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)

# opt=list()
# opt$citeSeq_data_path = "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsAnalyses/output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"

# Load pbmc data ----------------------------------------------------------

pbmc <- readRDS(opt$citeSeq_data_path)
pbmc$cell_type <- as.vector(pbmc$celltype.l2)

anno.table <- unique(data.frame("guiding_label" = pbmc$celltype.l1.5,
                                "fine_label" = pbmc$cell_type))
anno.table$guiding_label <- factor(anno.table$guiding_label,
                                   levels =c("B","CD4 T","CD8 T","dnT","gdT","MAIT","NK","ILC","Mono","DC","Eryth","Platelet","HSPC"))
anno.table <- unique(anno.table)
anno.table <- anno.table %>% arrange(factor(guiding_label),fine_label) 
anno.table

pbmc$cell_type <- factor(pbmc$cell_type,levels = anno.table$fine_label)


pbmc$time <- stringr::str_split_fixed(pbmc$orig.ident,pattern = "_",2)[,2]
pbmc$donor <- stringr::str_split_fixed(pbmc$orig.ident,pattern = "_",2)[,1]
DimPlot(pbmc,group.by = "celltype.l1",reduction = "wnn.umap")
pbmc <- FindMultiModalNeighbors(pbmc,reduction.list = list("pca","apca"),dims.list = list(c(1:40),c(1:50)),prune.SNN = 0)
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.uma", reduction.key = "wnnUMAP_",min.dist = 0.2,return.model = T) # go with min.dist=0.2 in fig 4
DimPlot(pbmc,reduction = "wnn.umap",label = T,group.by = "celltype.l2")

DefaultAssay(pbmc) <- "RNA"


pbmc <- FindClusters(pbmc,graph.name = "wsnn",algorithm = 3)
pbmc$metacell_cluster <- pbmc$wsnn_res.0.8

color_celltype <- colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442","#CCEBC5", "#D55E00", "#CC79A7",
  "#999999", "#AD7FA8", "#F781BF", "#A6D854", "#FFD92F", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#0072B2", "#CAB2D6", "#6A3D9A", "#B15928", "#BC80BD",
  "#FF7F00", "#B3B3B3", "#117733"
)
names(color_celltype) <- unique(pbmc$cell_type)
color_celltype_final <- color_celltype
color_celltype_final[["CD14 Mono"]] <- color_celltype[["CD8 Naive"]]
color_celltype_final[["CD8 Naive"]] <- color_celltype[["CD14 Mono"]]
# color_celltype_final[["CD4 TCM"]] <- color_celltype[["T/NK Prolif."]]
# color_celltype_final[["T/NK Prolif."]] <- color_celltype[["CD4 TCM"]]
color_celltype_final[["Platelet"]] <- color_celltype[["cDC1"]]
color_celltype_final[["cDC1"]] <- color_celltype[["Platelet"]]
color_celltype_final.2 <- color_celltype_final

# names(color_celltype_final.2) <- c(names(color_celltype_final.2)[1:28],"CD14 Mono/Platelet","CD14 Mono ifn")
pbmc$refined_celltype <- as.vector(pbmc$cell_type)
pbmc$refined_celltype[pbmc$wsnn_res.0.8 == 16] <- "CD14 Mono Ifn"
color_celltype_final.2[["CD14 Mono Ifn"]] <- "#b20000" 
DimPlot(pbmc,reduction = "wnn.umap",label = T,group.by = "refined_celltype",cols = color_celltype_final.2)


pbmc$metacell_cluster <- plyr::mapvalues(pbmc$wsnn_res.0.8,from = c(0,6,16),
                                         to = c("CD14 Mono","CD16 Mono", "CD14 Mono Ifn"))

pbmc$metacell_cluster <- factor(pbmc$metacell_cluster,levels =c("CD14 Mono","CD16 Mono", "CD14 Mono Ifn"))


# Select monocytes data ---------------------------------------------------

mono.pbmc <- pbmc[,pbmc$wsnn_res.0.8 %in% c(0,16,6)]
mono.pbmc <-  mono.pbmc[,mono.pbmc$metacell_cluster %in% c("CD14 Mono","CD16 Mono","CD14 Mono Ifn")]


# Run differential analysis using edgeR -----------------------------------

edgeR.obj <- Seurat2PB.custom(object = mono.pbmc, sample="orig.ident", 
                              assay = "RNA", cluster="metacell_cluster")
edgeR.obj$samples$cluster <- sapply(edgeR.obj$samples$cluster, function(i)gsub(x=i,pattern = " ",replacement = "_"))

summary(edgeR.obj$samples$lib.size)
table(edgeR.obj$samples$cluster)

edgeR.obj$samples$time <- sapply(edgeR.obj$samples$sample, function(i) unique(mono.pbmc$time[mono.pbmc$orig.ident == i]))
edgeR.obj$samples$donor <- sapply(edgeR.obj$samples$sample, function(i) unique(mono.pbmc$donor[mono.pbmc$orig.ident == i]))


keep.samples <- edgeR.obj$samples$lib.size > 4e4
table(keep.samples)
edgeR.obj <- edgeR.obj[, keep.samples]
table(edgeR.obj$samples$cluster)

keep.CPM <- rowSums(cpm(edgeR.obj) >= 10) >= (min(table(edgeR.obj$samples$cluster)* 0.7) )
table(keep.CPM)
keep.TotalCount <- (rowSums(edgeR.obj$counts) >= 10)
keep.genes <- keep.CPM & keep.TotalCount
table(keep.genes)

edgeR.obj <- edgeR.obj[keep.genes, , keep=FALSE]
edgeR.obj <- normLibSizes(edgeR.obj)

rna.genes <- rownames(edgeR.obj)


cluster <- factor(edgeR.obj$samples$cluster,levels =c("CD14_Mono","CD16_Mono", "CD14_Mono_Ifn"))
color.mds <- c("#56B4E9","#E69F00" ,"#009E73")
names(color.mds)<- levels(cluster)
plotMDS(edgeR.obj, pch=16, col=color.mds[cluster], main="MDS", top = 5000, gene.selection="common",) #, top = 5000, gene.selection="common"
legend("bottomleft", legend=levels(cluster), pch=16, cex=0.8)

time <- factor(edgeR.obj$samples$time)
color.mds <- c("#56B4E9","#E69F00" ,"#009E73")

names(color.mds)<- levels(time)
plotMDS(edgeR.obj, pch=16, col=color.mds[time], main="MDS", top = 5000, gene.selection="common",) #, top = 5000, gene.selection="common"
legend("bottomleft", legend=levels(time), pch=16, cex=0.8)


sample <- factor(edgeR.obj$samples$sample)
time <- factor(edgeR.obj$samples$time)
donor <- factor(edgeR.obj$samples$donor)


# design <- model.matrix(~ 0 + cluster*time + donor  ) 
design <- model.matrix(~ 0 + cluster:time + donor) 
# design <- model.matrix(~ 0 + time:cluster + donor) 
# design <- model.matrix(~ donor + cluster*time ) 
# design <- model.matrix(~ donor + cluster ) 
# design <- model.matrix(~ donor + cluster:time )
# design <- model.matrix(~ 0 + time:cluster + donor)
# design <- model.matrix(~ 0 + donor + time:cluster )
# design <- model.matrix(~ 0 + cluster + time + donor) ??
# design <- model.matrix(~ 0 + cluster * time + donor) ??

# colnames(design)[1] <- "Int"
head(design)

edgeR.obj <- estimateDisp(edgeR.obj, design, robust=TRUE)
edgeR.obj$common.dispersion
gc()
plotBCV(edgeR.obj)

fit <- glmQLFit(edgeR.obj, design, robust=TRUE)
gc()
plotQLDisp(fit)


# Define contrasts for each cluster
colnames(design) <- make.names(colnames(design))

contrasts <- makeContrasts(
  CD14_Mono_vs_All = clusterCD14_Mono - (clusterCD16_Mono + clusterCD14_Mono_Ifn) / 2,
  CD16_Mono_vs_All = clusterCD16_Mono - (clusterCD14_Mono + clusterCD14_Mono_Ifn) / 2,
  CD14_Mono_Ifn_vs_All = clusterCD14_Mono_Ifn - (clusterCD14_Mono + clusterCD16_Mono) / 2,
  levels = design
)

# Test for Cluster1 markers
results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_All"])


# Extract significant markers for each cluster
res.all.genes <- topTags(results, n = Inf, sort.by = "PValue")$table
write.table(res.all.genes, file = paste0(opt$output_path, "/edgeR_res_all.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
res <- res.all.genes[(res.all.genes$FDR < 0.05 & res.all.genes$logFC >0.25)&res.all.genes$logCPM>4,]
res


# Pairwise comparisons ----------------------------------------------------


contrasts <- makeContrasts(
  CD14_Mono_Ifn_vs_CD14_Mono = clusterCD14_Mono_Ifn - clusterCD14_Mono,
  CD14_Mono_Ifn_vs_CD16_Mono = clusterCD14_Mono_Ifn - clusterCD16_Mono,
  levels = design
)

results_vs_CD14 <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_CD14_Mono"])
results_vs_CD16 <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_CD16_Mono"])


# Extract significant markers for each cluster
res_vs_CD14 <- topTags(results_vs_CD14, n = Inf, sort.by = "PValue")$table
res_vs_CD16 <- topTags(results_vs_CD16, n = Inf, sort.by = "PValue")$table
colnames(res_vs_CD14) <- paste0(colnames(res_vs_CD14) ,"_vs_CD14")
colnames(res_vs_CD16) <- paste0(colnames(res_vs_CD16) ,"_vs_CD16")
common.genes <- rownames(res_vs_CD14)[rownames(res_vs_CD14) %in% rownames(res_vs_CD16)]

res_pairwise  <- cbind(res_vs_CD14[common.genes,],res_vs_CD16[common.genes,])
res_pairwise <- res_pairwise[res_pairwise$logCPM_vs_CD14*res_pairwise$logCPM_vs_CD16 >0,]
res_pairwise$max_FDR <- apply(res_pairwise[,c("FDR_vs_CD14","FDR_vs_CD16")],MARGIN =1, FUN = max)
res_pairwise <- res_pairwise[order(res_pairwise$max_FDR,decreasing = F),]
res_pairwise <- res_pairwise[(res_pairwise$max_FDR< 0.05 & res_pairwise$logCPM_vs_CD14>4)&res_pairwise$logCPM_vs_CD16>4,]
res_pairwise
write.table(res_pairwise, file = paste0(opt$output_path, "/edgeR_res_pairwise.txt"), quote = F, row.names = F, col.names = T, sep = "\t")