
library(edgeR)
library(getopt)
library(dplyr)
library(Seurat)
print(getwd())
source("figures/fun_r/seurat2PB.R")

spec <- matrix(c(
  'citeSeq_data_path',  'c', 1, "character",
  'output_path',  'o', 1, "character",
  'help',   'h', 0, "logical",
  'diff_assay',  'a', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)

if(is.null(opt$diff_assay)){
  opt$diff_assay <- "RNA"
}
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
print(mono.pbmc)
print(opt$diff_assay)
edgeR.obj <- Seurat2PB.custom(object = mono.pbmc, sample="orig.ident", 
                              assay = opt$diff_assay, cluster="metacell_cluster")
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
plotMDS(edgeR.obj, pch=16, col=color.mds[cluster], main="MDS", top = 5000, gene.selection="common",dim.plot = 1:2) #, top = 5000, gene.selection="common"
legend("bottomleft", legend=levels(cluster), pch=16, cex=0.8)

plotMDS(edgeR.obj, pch=16, col=color.mds[cluster], main="MDS", top = 5000, gene.selection="common",dim.plot = c(1,3)) #, top = 5000, gene.selection="common"
legend("bottomleft", legend=levels(cluster), pch=16, cex=0.8)

time <- factor(edgeR.obj$samples$time)
color.mds <- c("#56B4E9","#E69F00" ,"#009E73")

names(color.mds)<- levels(time)
plotMDS(edgeR.obj, pch=16, col=color.mds[time], main="MDS", top = 5000, gene.selection="common",) #, top = 5000, gene.selection="common"
legend("bottomleft", legend=levels(time), pch=16, cex=0.8)

# color.mds <-color_celltype
# names(color.mds)<- levels(donor)
# plotMDS(edgeR.obj, pch=16, col=color.mds[donor], main="MDS", top = 5000, gene.selection="common",) #, top = 5000, gene.selection="common"
# legend("bottomleft", legend=levels(donor), pch=16, cex=0.8)

sample <- factor(edgeR.obj$samples$sample)
time <- factor(edgeR.obj$samples$time)
donor <- factor(edgeR.obj$samples$donor)

Group <- factor(paste(cluster,time,sep="."), levels = unique(paste(cluster,time,sep=".")))


###########################################################################
###########################################################################
###                                                                     ###
###                            ALL MONOCYTES                            ###
###                                                                     ###
###########################################################################
###########################################################################

# design <- model.matrix(~ 0 + cluster*time + donor  )
design <- model.matrix(~ donor + Group)
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
  # differences of each cell type vs others at T0
  CD14_Mono_vs_All_T0     = -(GroupCD16_Mono.0 + GroupCD14_Mono_Ifn.0)/2,
  CD16_Mono_vs_All_T0     =  GroupCD16_Mono.0 - (GroupCD14_Mono_Ifn.0)/2,
  CD14_Mono_Ifn_vs_All_T0 =  GroupCD14_Mono_Ifn.0 - (GroupCD16_Mono.0)/2,
  # differences of each cell type vs others at T2
  CD14_Mono_vs_All_T2     =  GroupCD14_Mono.2     - (GroupCD16_Mono.2     + GroupCD14_Mono_Ifn.2)/2,
  CD16_Mono_vs_All_T2     =  GroupCD16_Mono.2     - (GroupCD14_Mono.2     + GroupCD14_Mono_Ifn.2)/2,
  CD14_Mono_Ifn_vs_All_T2 =  GroupCD14_Mono_Ifn.2 - (GroupCD14_Mono.2     + GroupCD16_Mono.2)/2,
  # differences of each cell type vs others at T7
  CD14_Mono_vs_All_T7     =  GroupCD14_Mono.7     - (GroupCD16_Mono.7     + GroupCD14_Mono_Ifn.7)/2,
  CD16_Mono_vs_All_T7     =  GroupCD16_Mono.7     - (GroupCD14_Mono.7     + GroupCD14_Mono_Ifn.7)/2,
  CD14_Mono_Ifn_vs_All_T7 =  GroupCD14_Mono_Ifn.7 - (GroupCD14_Mono.7     + GroupCD16_Mono.7)/2,
  levels = design
)

# Extract significant markers
results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_All_T0"])
res.all.genes.t0 <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes.t0$timepoint <- "t0"

res.t0 <- res.all.genes.t0[(res.all.genes.t0$FDR < 0.05 & abs(res.all.genes.t0$logFC) >0.25) & res.all.genes.t0$logCPM>4,]
res.t0 <- res.t0[order(res.t0$logFC, decreasing = T), ]
dim(res.t0)
write.table(res.all.genes.t0, file = paste0(opt$output_path, "/edgeR_res_all_t0.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_All_T2"])
res.all.genes.t2 <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes.t2$timepoint <- "t2"

res.t2 <- res.all.genes.t2[(res.all.genes.t2$FDR < 0.05 & abs(res.all.genes.t2$logFC) >0.25) & res.all.genes.t2$logCPM>4,]
res.t2 <- res.t2[order(res.t2$logFC, decreasing = T), ]
dim(res.t2)
write.table(res.all.genes.t2, file = paste0(opt$output_path, "/edgeR_res_all_t2.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_All_T7"])
res.all.genes.t7 <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes.t7$timepoint <- "t7"

res.t7 <- res.all.genes.t7[(res.all.genes.t7$FDR < 0.05 & abs(res.all.genes.t7$logFC) >0.25) & res.all.genes.t7$logCPM>4,]
res.t7 <- res.t7[order(res.t7$logFC, decreasing = T), ]
dim(res.t7)
write.table(res.all.genes.t7, file = paste0(opt$output_path, "/edgeR_res_all_t7.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# save most differential genes in all timepoints in a file
combined.res <- rbind(res.t0, res.t2, res.t7)
combined.res$annotation <- "other"
# save timepoint specific genes
t0.specific <- res.t0$gene[!res.t0$gene %in% c(res.t2$gene, res.t7$gene)]
length(t0.specific)
t2.specific <- res.t2$gene[!res.t2$gene %in% c(res.t0$gene, res.t7$gene)]
length(t2.specific)
t7.specific <- res.t7$gene[!res.t7$gene %in% c(res.t2$gene, res.t0$gene)]
length(t7.specific)

# identify the DEGs in common between the 3 timepoints
common_DEGs <- Reduce(intersect, list(
  res.t0$gene,
  res.t2$gene,
  res.t7$gene
))
length(common_DEGs)
combined.res$annotation[combined.res$gene %in% t0.specific] <- "t0_specific"
combined.res$annotation[combined.res$gene %in% t2.specific] <- "t2_specific"
combined.res$annotation[combined.res$gene %in% t7.specific] <- "t7_specific"
combined.res$annotation[combined.res$gene %in% common_DEGs] <- "common"

write.table(combined.res, file = paste0(opt$output_path, "/edgeR_res_all_summary.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# Get genes differentially expressed between CD14_Ifn and the other monocytes at any timepoint
results <- glmQLFTest(fit, contrast = contrasts[,c("CD14_Mono_Ifn_vs_All_T0", "CD14_Mono_Ifn_vs_All_T2", "CD14_Mono_Ifn_vs_All_T7")])
res.all.genes <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes$minFC <- matrixStats::rowMins(abs(as.matrix(res.all.genes[, c("logFC.CD14_Mono_Ifn_vs_All_T0", "logFC.CD14_Mono_Ifn_vs_All_T2", "logFC.CD14_Mono_Ifn_vs_All_T7")])))
res.all.genes$meanFC <- rowMeans(as.matrix(res.all.genes[, c("logFC.CD14_Mono_Ifn_vs_All_T0", "logFC.CD14_Mono_Ifn_vs_All_T2", "logFC.CD14_Mono_Ifn_vs_All_T7")]))

res.allT <- res.all.genes[(res.all.genes$FDR < 0.05 & abs(res.all.genes$minFC) >0.25)&res.all.genes$logCPM>4,]
res.allT <- res.allT[order(res.allT$meanFC, decreasing = T), ]
length(res.allT$gene)

# check that most of genes differentially expressed in any of the 3 timepoints are actually aprt of the DEGs found at a particular time point
summary(res.allT$gene %in% c(res.t0$gene, res.t2$gene, res.t7$gene))

write.table(res.all.genes, file = paste0(opt$output_path, "/edgeR_res_all.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# save edgeR data 
lcpm <- cpm(edgeR.obj, log=TRUE, normalized.lib.sizes = T)
annot <- data.frame(cluster=cluster,
                    timepoint = paste0("timepoint_", time))
rownames(annot) <- colnames(edgeR.obj)
save(lcpm, annot, file = paste0(opt$output_path, "edgeR_data.rdata"))

# Pairwise comparisons ----------------------------------------------------

contrasts <- makeContrasts(
  CD14_Mono_Ifn_vs_CD14_Mono = GroupCD14_Mono_Ifn.0,
  CD14_Mono_Ifn_vs_CD16_Mono = GroupCD14_Mono_Ifn.0 - GroupCD16_Mono.0,
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


###########################################################################
###########################################################################
###                                                                     ###
###                         ONLY CD14 MONOCYTES                         ###
###                                                                     ###
###########################################################################
###########################################################################

edgeR.obj <- edgeR.obj[, cluster != "CD16_Mono"]

cluster <- factor(edgeR.obj$samples$cluster,levels =c("CD14_Mono", "CD14_Mono_Ifn"))
sample <- factor(edgeR.obj$samples$sample)
time <- factor(edgeR.obj$samples$time)
donor <- factor(edgeR.obj$samples$donor)

Group <- factor(paste(cluster,time,sep="."), levels = unique(paste(cluster,time,sep=".")))

design <- model.matrix(~ donor + Group)
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
  # differences of each cell type vs others at T0
  CD14_Mono_Ifn_vs_CD14_Mono_T0 =  GroupCD14_Mono_Ifn.0,
  # differences of each cell type vs others at T2
  CD14_Mono_Ifn_vs_CD14_Mono_T2 =  GroupCD14_Mono_Ifn.2 - GroupCD14_Mono.2,
  # differences of each cell type vs others at T7
  CD14_Mono_Ifn_vs_CD14_Mono_T7 = GroupCD14_Mono_Ifn.7 - GroupCD14_Mono.7,
  levels = design
)

# Extract significant markers
results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_CD14_Mono_T0"])
res.all.genes.t0 <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes.t0$timepoint <- "t0"

res.t0 <- res.all.genes.t0[(res.all.genes.t0$FDR < 0.05 & abs(res.all.genes.t0$logFC) >0.25) & res.all.genes.t0$logCPM>4,]
res.t0 <- res.t0[order(res.t0$logFC, decreasing = T), ]
dim(res.t0)
write.table(res.all.genes.t0, file = paste0(opt$output_path, "/edgeR_res_CD14_t0.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_CD14_Mono_T2"])
res.all.genes.t2 <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes.t2$timepoint <- "t2"

res.t2 <- res.all.genes.t2[(res.all.genes.t2$FDR < 0.05 & abs(res.all.genes.t2$logFC) >0.25) & res.all.genes.t2$logCPM>4,]
res.t2 <- res.t2[order(res.t2$logFC, decreasing = T), ]
dim(res.t2)
write.table(res.all.genes.t2, file = paste0(opt$output_path, "/edgeR_res_CD14_t2.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

results <- glmQLFTest(fit, contrast = contrasts[,"CD14_Mono_Ifn_vs_CD14_Mono_T7"])
res.all.genes.t7 <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes.t7$timepoint <- "t7"

res.t7 <- res.all.genes.t7[(res.all.genes.t7$FDR < 0.05 & abs(res.all.genes.t7$logFC) >0.25) & res.all.genes.t7$logCPM>4,]
res.t7 <- res.t7[order(res.t7$logFC, decreasing = T), ]
dim(res.t7)
write.table(res.all.genes.t7, file = paste0(opt$output_path, "/edgeR_res_CD14_t7.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# save most differential genes in all timepoints in a file
combined.res <- rbind(res.t0, res.t2, res.t7)
combined.res$annotation <- "other"
# save timepoint specific genes
t0.specific <- res.t0$gene[!res.t0$gene %in% c(res.t2$gene, res.t7$gene)]
length(t0.specific)
t2.specific <- res.t2$gene[!res.t2$gene %in% c(res.t0$gene, res.t7$gene)]
length(t2.specific)
t7.specific <- res.t7$gene[!res.t7$gene %in% c(res.t2$gene, res.t0$gene)]
length(t7.specific)

# identify the DEGs in common between the 3 timepoints
common_DEGs <- Reduce(intersect, list(
  res.t0$gene,
  res.t2$gene,
  res.t7$gene
))
length(common_DEGs)
combined.res$annotation[combined.res$gene %in% t0.specific] <- "t0_specific"
combined.res$annotation[combined.res$gene %in% t2.specific] <- "t2_specific"
combined.res$annotation[combined.res$gene %in% t7.specific] <- "t7_specific"
combined.res$annotation[combined.res$gene %in% common_DEGs] <- "common"

write.table(combined.res, file = paste0(opt$output_path, "/edgeR_res_CD14_summary.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# Get genes differentially expressed between CD14_Ifn and the other monocytes at any timepoint
results <- glmQLFTest(fit, contrast = contrasts[,c("CD14_Mono_Ifn_vs_CD14_Mono_T0", "CD14_Mono_Ifn_vs_CD14_Mono_T2", "CD14_Mono_Ifn_vs_CD14_Mono_T7")])
res.all.genes <- topTags(results, n = Inf, sort.by = "PValue")$table
res.all.genes$minFC <- matrixStats::rowMins(abs(as.matrix(res.all.genes[, c("logFC.CD14_Mono_Ifn_vs_CD14_Mono_T0", "logFC.CD14_Mono_Ifn_vs_CD14_Mono_T2", "logFC.CD14_Mono_Ifn_vs_CD14_Mono_T7")])))
res.all.genes$meanFC <- rowMeans(as.matrix(res.all.genes[, c("logFC.CD14_Mono_Ifn_vs_CD14_Mono_T0", "logFC.CD14_Mono_Ifn_vs_CD14_Mono_T2", "logFC.CD14_Mono_Ifn_vs_CD14_Mono_T7")]))

res.allT <- res.all.genes[(res.all.genes$FDR < 0.05 & abs(res.all.genes$minFC) >0.25)&res.all.genes$logCPM>4,]
res.allT <- res.allT[order(res.allT$meanFC, decreasing = T), ]
length(res.allT$gene)

# check that most of genes differentially expressed in any of the 3 timepoints are actually aprt of the DEGs found at a particular time point
summary(res.allT$gene %in% c(res.t0$gene, res.t2$gene, res.t7$gene))

write.table(res.all.genes, file = paste0(opt$output_path, "/edgeR_res_CD14.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# save edgeR data 
lcpm <- cpm(edgeR.obj, log=TRUE, normalized.lib.sizes = T)
annot <- data.frame(cluster=cluster,
                    timepoint = paste0("timepoint_", time))
rownames(annot) <- colnames(edgeR.obj)
save(lcpm, annot, file = paste0(opt$output_path, "edgeR_data_CD14.rdata"))


