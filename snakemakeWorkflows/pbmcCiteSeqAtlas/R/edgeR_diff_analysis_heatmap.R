library(edgeR)
library(getopt)
library(dplyr)
library(Seurat)

opt=list()
opt$citeSeq_data_path = "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsAnalyses/output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
opt$output_path = "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsAnalyses/output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/"

###########################################################################
###########################################################################
###                                                                     ###
###                            ALL MONOCYTES                            ###
###                                                                     ###
###########################################################################
###########################################################################

res.all.genes.t0 <- read.table(paste0(opt$output_path, "/edgeR_res_all_t0.txt"), header = T, sep = "\t")
res.t0 <- res.all.genes.t0[(res.all.genes.t0$FDR < 0.05 & res.all.genes.t0$logFC >0.25) & res.all.genes.t0$logCPM>4,]
res.t0 <- res.t0[order(res.t0$logFC, decreasing = T), ]
dim(res.t0)

res.all.genes.t2 <- read.table(paste0(opt$output_path, "/edgeR_res_all_t2.txt"), header = T, sep = "\t")
res.t2 <- res.all.genes.t2[(res.all.genes.t2$FDR < 0.05 & res.all.genes.t2$logFC >0.25) & res.all.genes.t2$logCPM>4,]
res.t2 <- res.t2[order(res.t2$logFC, decreasing = T), ]
dim(res.t2)

res.all.genes.t7 <- read.table(paste0(opt$output_path, "/edgeR_res_all_t7.txt"), header = T, sep = "\t")
res.t7 <- res.all.genes.t7[(res.all.genes.t7$FDR < 0.05 & res.all.genes.t7$logFC >0.25) & res.all.genes.t7$logCPM>4,]
res.t7 <- res.t7[order(res.t7$logFC, decreasing = T), ]
dim(res.t7)


res.pairwise <- read.table(paste0(opt$output_path, "/edgeR_res_pairwise.txt"), header = T, sep = "\t")
res.pairwise <- res.pairwise[res.pairwise$FDR_vs_CD14 < 0.05 & res.pairwise$FDR_vs_CD16 < 0.05 &
                                res.pairwise$logFC_vs_CD14 >0.25 & res.pairwise$logFC_vs_CD16 >0.25 & res.pairwise$logCPM_vs_CD14>4,]
res.pairwise <- res.pairwise[order(res.pairwise$max_FDR, decreasing = T), ]
dim(res.pairwise)

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
common_FC <- sapply(common_DEGs, function(i) max(combined.res$logFC[combined.res$gene == i]))
common_DEGs <- common_DEGs[order(common_FC, decreasing = T)]

# Visualization -----------------------------------------------------------
load(paste0(opt$output_path, "edgeR_data.rdata"))
ann_colors <- list(cluster=c("CD14_Mono" = "#56B4E9",  "CD16_Mono" = "#E69F00" , "CD14_Mono_Ifn" = "#009E73"),
                   timepoint = c("timepoint_0" = '#E0D4CA', "timepoint_2" = '#476D87', "timepoint_7" = '#B53E2B'),
                   group = c("t0_specific"='#D35FB7',"t2_specific" = '#0072B2', "t7_specific" = '#00BFC4', "common" = '#FF6F00'))

ordered_cols <- annot %>%
  arrange(timepoint, cluster) 
gaps_c <- cumsum(table(annot$timepoint))

pheatmap::pheatmap(lcpm[res.t0$gene[1:100], rownames(ordered_cols)], 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   scale="row",
                   cluster_cols=F, cluster_rows = T, 
                   border_color="NA", fontsize_row=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors,
                   gaps_col = gaps_c)


pheatmap::pheatmap(lcpm[res.pairwise$gene_vs_CD14, rownames(ordered_cols)], 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   scale="row",
                   cluster_cols=F, cluster_rows = T, 
                   border_color="NA", fontsize_row=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors,
                   gaps_col = gaps_c)


ntop = 20
topMarkers <- c(t0.specific[1:ntop], t2.specific[1:ntop], t7.specific[1:ntop], common_DEGs[1:ntop])
gene_group = data.frame(row.names = topMarkers, 
                        group = c(rep("t0_specific", ntop), 
                                  rep("t2_specific", ntop), 
                                  rep("t7_specific", ntop), 
                                  rep("common", ntop)))
gene_group$group <- factor(gene_group$group, levels = c("t0_specific", "t2_specific", "t7_specific", "common"))

gaps <- cumsum(table(gene_group$group))
gaps_c <- cumsum(table(annot$timepoint))


pheatmap::pheatmap(lcpm[topMarkers, rownames(ordered_cols)], breaks=seq(-4,4,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   scale="row",
                   cluster_cols=F, cluster_rows = F, 
                   border_color="NA", fontsize_row=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors,
                   annotation_row = gene_group,
                   gaps_row = gaps,
                   gaps_col = gaps_c)


pheatmap::pheatmap(lcpm[t0.specific[1:100], rownames(ordered_cols)], breaks=seq(-4,4,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   scale="row",
                   cluster_cols=F, cluster_rows = F, 
                   border_color="NA", fontsize_row=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors,
                   annotation_row = gene_group,
                   gaps_row = gaps,
                   gaps_col = gaps_c)



###########################################################################
###########################################################################
###                                                                     ###
###                         ONLY CD14 MONOCYTES                         ###
###                                                                     ###
###########################################################################
###########################################################################

res.all.genes.t0 <- read.table(paste0(opt$output_path, "/edgeR_res_CD14_t0.txt"), header = T, sep = "\t")
res.t0 <- res.all.genes.t0[(res.all.genes.t0$FDR < 0.05 & res.all.genes.t0$logFC >0.25) & res.all.genes.t0$logCPM>4,]
res.t0 <- res.t0[order(res.t0$logFC, decreasing = T), ]
dim(res.t0)

res.all.genes.t2 <- read.table(paste0(opt$output_path, "/edgeR_res_CD14_t2.txt"), header = T, sep = "\t")
res.t2 <- res.all.genes.t2[(res.all.genes.t2$FDR < 0.05 & res.all.genes.t2$logFC >0.25) & res.all.genes.t2$logCPM>4,]
res.t2 <- res.t2[order(res.t2$logFC, decreasing = T), ]
dim(res.t2)

res.all.genes.t7 <- read.table(paste0(opt$output_path, "/edgeR_res_CD14_t7.txt"), header = T, sep = "\t")
res.t7 <- res.all.genes.t7[(res.all.genes.t7$FDR < 0.05 & res.all.genes.t7$logFC >0.25) & res.all.genes.t7$logCPM>4,]
res.t7 <- res.t7[order(res.t7$logFC, decreasing = T), ]
dim(res.t7)


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
common_FC <- sapply(common_DEGs, function(i) max(combined.res$logFC[combined.res$gene == i]))
common_DEGs <- common_DEGs[order(common_FC, decreasing = T)]

# Visualization -----------------------------------------------------------
load(paste0(opt$output_path, "edgeR_data_CD14.rdata"))
ann_colors <- list(cluster=c("CD14_Mono" = "#56B4E9",  "CD16_Mono" = "#E69F00" , "CD14_Mono_Ifn" = "#009E73"),
                   timepoint = c("timepoint_0" = '#E0D4CA', "timepoint_2" = '#476D87', "timepoint_7" = '#B53E2B'),
                   group = c("t0_specific"='#D35FB7',"t2_specific" = '#0072B2', "t7_specific" = '#00BFC4', "common" = '#FF6F00'))

ordered_cols <- annot %>%
  arrange(timepoint, cluster) 
gaps_c <- cumsum(table(annot$timepoint))

pheatmap::pheatmap(lcpm[res.t0$gene[1:100], rownames(ordered_cols)], 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   scale="row",
                   cluster_cols=F, cluster_rows = T, 
                   border_color="NA", fontsize_row=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors,
                   gaps_col = gaps_c)


ntop = 20
topMarkers <- c(t0.specific[1:min(length(t0.specific),ntop)], t2.specific[1:min(length(t2.specific),ntop)], 
                t7.specific[1:min(length(t7.specific),ntop)], common_DEGs[1:min(length(common_DEGs),ntop)])
gene_group = data.frame(row.names = topMarkers, 
                        group = c(rep("t0_specific", min(length(t0.specific),ntop)), 
                                  rep("t2_specific", min(length(t2.specific),ntop)), 
                                  rep("t7_specific", min(length(t7.specific),ntop)), 
                                  rep("common", min(length(common_DEGs),ntop))))
gene_group$group <- factor(gene_group$group, levels = c("t0_specific", "t2_specific", "t7_specific", "common"))

gaps <- cumsum(table(gene_group$group))
gaps_c <- cumsum(table(annot$timepoint))


pheatmap::pheatmap(lcpm[topMarkers, rownames(ordered_cols)], breaks=seq(-4,4,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   scale="row",
                   cluster_cols=F, cluster_rows = F, 
                   border_color="NA", fontsize_row=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors,
                   annotation_row = gene_group,
                   gaps_row = gaps,
                   gaps_col = gaps_c)


# entrez_ids <- AnnotationDbi::mapIds(
#   org.Hs.eg.db::org.Hs.eg.db,
#   keys = common,             # your gene symbols
#   column = "ENTREZID",      # what you want
#   keytype = "SYMBOL",       # what you have
#   multiVals = "first"       # if multiple matches, take first
# )
# go <- goana(entrez_ids, species="Hs")
# go <- go[go$Ont == "BP" & go$P.DE < 0.05,]

