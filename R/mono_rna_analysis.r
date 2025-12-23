# Libraries ----------------------------------------------------------------

library(edgeR)
library(tximport)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(GSVA)
library(ggpubr)
library(msigdbr)
library(patchwork)
library(ComplexHeatmap)
library(readxl)
library(rtracklayer)

# Parameters --------------------------------------------------------------

# output_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/downstream_analyses/"
# dir.create(output_path)
# dir.create(paste0(output_path, "/figures/"))
# input_data_salmon <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/preprocessing_results/star_salmon/"
# sample_sheet <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/samplesheet.csv"
# sinature_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsHTAN/data/"
# citeseq_mono_markers <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsAnalyses/output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_CD14_t0.txt"
# htan_mac_markers <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsHTAN/output/whole_atlas/metacells_sup/g10_final/marker_analysis/Mono_Macro/RNA_diff_res.csv"
# gtf_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/Homo_sapiens.GRCh38.115.gtf.gz"

output_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/downstream_analyses/"
dir.create(output_path)
dir.create(paste0(output_path, "/figures/"))
input_data_salmon <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/preprocessing_results/star_salmon/"
sample_sheet <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/samplesheet.csv"
sinature_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsHTAN/data/"
citeseq_mono_markers <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsAnalyses/output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_CD14_t0.txt"
htan_mac_markers <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/SuperCellMultiomicsHTAN/output/whole_atlas/metacells_sup/g10_final/marker_analysis/Mono_Macro/RNA_diff_res.csv"
gtf_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/preprint_final/manuscript_results/RNA_data/Homo_sapiens.GRCh38.115.gtf.gz"


# Load metadata ---------------------------------------------------------

metadata <- read.csv2(sample_sheet, sep = ",")
metadata$donor <- gsub("_.*", "", metadata$sample)
metadata$group <- gsub(".*_", "", metadata$sample)
head(metadata)
table(metadata$group)
metadata$group <- factor(metadata$group, levels = c("DN", "CD169", "LY6E", "DP"))

# Load salmon quantification ----------------------------------------------

gtf <- import(gtf_path)
head(gtf)
protein_coding_genes <- gtf %>%
  as.data.frame() %>%
  filter(type %in% c("gene", 'transcript'), gene_biotype == "protein_coding")
# summary(protein_coding_genes$transcript_id %in% tx2gene$transcript_ID)

# load transcript annotation 
tx2gene <- read_table(paste0(input_data_salmon, "/tx2gene.tsv"), col_names = c("transcript_ID","gene_ID","gene_symbol"))
tx2gene %>% head()

gene_matching <- as.data.frame(unique(tx2gene[, 2:3]))
rownames(gene_matching) <- gene_matching$gene_ID
dim(gene_matching)

files <- file.path(input_data_salmon, metadata$sample, "quant.sf")
names(files) <- metadata$sample
files

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
head(txi)

# round the pseudocounts 
counts_salmon <- txi$counts %>% 
  round() %>% 
  data.frame()
dim(counts_salmon)

counts_with_names <- cbind(counts_salmon, gene_matching[rownames(counts_salmon), ])
counts_salmon <- counts_with_names %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
counts_salmon <- as.data.frame(counts_salmon)
rownames(counts_salmon) <- counts_salmon$gene_symbol
counts_salmon <- counts_salmon[, -1]
dim(counts_salmon)

tpm_salmon <- txi$abundance
rownames(metadata) <- metadata$sample

tpm_with_names <- cbind(tpm_salmon, gene_matching[rownames(tpm_salmon), ])
tpm_salmon <- tpm_with_names %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
tpm_salmon <- as.data.frame(tpm_salmon)
rownames(tpm_salmon) <- tpm_salmon$gene_symbol
tpm_salmon <- tpm_salmon[, -1]
dim(tpm_salmon)

# keep genes with minimal level of expression
keep <- rowSums(tpm_salmon > 1) >= 2
summary(keep)
tpm_salmon_filtered <- tpm_salmon[keep, ]


# Downstream analyses -----------------------------------------------------
counts <- counts_salmon
pdata <- counts %>% 
  gather(key = Sample, value = Count)
head(pdata)

ggplot(pdata) +
  geom_density(aes(x = Count, color = Sample)) +
  facet_wrap(~ Sample)+
  xlab("Raw expression counts") +
  ylab("Number of genes")



# Run PCA -----------------------------------------------------------------
tpm <- tpm_salmon_filtered

norm_data=log(tpm+1)

# Keep most variable peaks 
rvdm = apply(norm_data,1,var)
selectdm = order(rvdm, decreasing = TRUE)[1:500]
norm_data_subset = t(norm_data[selectdm,])
norm_data_subset[1:5,1:5]
dim(norm_data_subset)

# Run PCA and UMAP 
pca = prcomp(norm_data_subset, scale. = T)
print(summary(pca))
plot(pca)
pcaData = as.data.frame(pca$x)
loadings = pca$rotation
rownames(pcaData)=rownames(norm_data_subset)

to_plot=merge(pcaData, metadata, by.x="row.names",by.y="sample")
pdf(paste0(output_path,"figures/pca.pdf"), w=7, h=5)
ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape = donor)) +
  geom_point(size=2) + scale_shape_manual(values = c(16,17,18,1,2,3,7,9)) +
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Group",override.aes = list(size=4)))

ggplot(to_plot, aes(x=PC1, y=PC3, color=group,shape = donor)) +
  geom_point(size=2) + scale_shape_manual(values = c(16,17,18,1,2,3,7,9)) +
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Group",override.aes = list(size=4)))
dev.off()

pdf(paste0(output_path, "figures/marker_levels.pdf"))
for(gene in c("LY6E", "CD14","FCGR3A", "SIGLEC1", "CXCL10")){
  pData <- as.data.frame(t(tpm))
  pData$group <- factor(metadata$group, levels = c("DN", "LY6E", "CD169", "DP"))
  pData$sample <- metadata$sample
  pData$donor <- metadata$donor

  
  print(ggplot(pData, aes(x = group, y = .data[[gene]])) +
          geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0)) +
          geom_line(aes(group = donor, color = donor), linewidth = 0.7, alpha = 0.7) +
          stat_summary(fun = median, geom = "crossbar",
                       width = 0.4, fatten = 2, color = "black", linewidth = 0.5) +
          stat_compare_means(comparisons =combn(levels(factor(pData$group)), 2, simplify = FALSE),
                             method = "wilcox.test",
                             paired = T,       
                             label = "p.format")  +
          labs(title = gene, x = "Group", y = "Expression (TPM)") +
          theme_minimal() +
          theme())

}
dev.off()



# Differential analysis on the signatures ---------------------------------

split.vector <- factor(metadata[, "group"], 
                       levels = c("DN",  "LY6E","CD169", "DP"))
conditions <- metadata$group
comp_list = list(c("DN", "CD169"),c("DN", "LY6E"), c("DN", "DP"), c("CD169", "LY6E"),
                 c("CD169", "DP"),c("LY6E", "DP"))

descriminative_colors <- c("DN"="#ccebc5ff", "CD169" = "#6a1b9aff", "LY6E" = "#009688ff",
                           "DP" = "#b20000ff")


donor <- metadata$donor
groups <- unique(conditions)

design <- model.matrix(~donor+conditions)

y <- DGEList(counts=counts_salmon, samples = metadata)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

contrasts <- makeContrasts(
  CD169_vs_DN     = conditionsCD169,
  LY6E_vs_DN     =  conditionsLY6E,
  DP_vs_DN =  conditionsDP,
  LY6E_vs_CD169     =  conditionsLY6E - conditionsCD169,
  DP_vs_CD169     =  conditionsDP - conditionsCD169,
  DP_vs_LY6E =  conditionsDP - conditionsLY6E,
  levels = design
)


# Extract significant markers
summary_df_degs <- vector()
for(c in colnames(contrasts)){
  qlf <- glmQLFTest(fit, contrast = contrasts[,c])
  res.all <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  res.all$group1 <- strsplit(c,"_vs_")[[1]][1]
  res.all$group2 <- strsplit(c,"_vs_")[[1]][2]
  res.all$gene_symbol <- rownames(res.all)
  summary_df_degs <- rbind(summary_df_degs, res.all)
}

summary_df_degs_signif <- summary_df_degs[summary_df_degs$FDR <= 0.05, ]
DEGs <- unique(summary_df_degs_signif$gene_symbol)
# DEGs2 <- unique(summary_df_degs_signif$gene_symbol[summary_df_degs_signif$group1 == "DP" & summary_df_degs_signif$group2 == "DN"])
# summary(DEGs %in% DEGs2)
# DEGs[!DEGs %in% DEGs2]

DEGs <- DEGs[DEGs %in% protein_coding_genes$gene_name]

colours <- list('Condition' = descriminative_colors)
colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
                                            which = 'col',
                                            col = colours,
                                            annotation_legend_param = list(
                                              Condition = list(title = "Cell type")
                                            ),show_annotation_name = FALSE)



get_pval <- function(group1, group2, deg_res, genes){
  deg_res <- deg_res[which(deg_res$group1 == group1 & deg_res$group2 == group2), ]
  rownames(deg_res) <- deg_res$gene_symbol
  pval <- deg_res[genes, "FDR"]
  names(pval) = genes
  
  is_sig = pval < 0.05
  pch = rep("*", length(pval))
  pch[!is_sig] = NA#"ns"
  pch[pval<0.01] = "**"
  pch[pval<0.001] = "***"
  
  return(list(pval, pch))
}

pvalue_col_fun = circlize::colorRamp2(c(0, 2, 3), viridis::viridis(3)) 
rowAnn = rowAnnotation(
  `LY6E-,CD169-` = anno_simple(-log10(get_pval(group1 = "DP", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[1]]), 
                       col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval(group1 = "DP", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[2]]),
  `LY6E+,CD169-` = anno_simple(-log10(get_pval(group1 = "DP", group2 = "LY6E", deg_res = summary_df_degs, genes = DEGs)[[1]]), 
                           col = pvalue_col_fun,gp = gpar(col = "black"), pch = get_pval(group1 = "DP", group2 = "LY6E", deg_res = summary_df_degs, genes = DEGs)[[2]]),
  `LY6E-,CD169+` = anno_simple(-log10(get_pval(group1 = "DP", group2 = "CD169", deg_res = summary_df_degs, genes = DEGs)[[1]]), 
                         col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval(group1 = "DP", group2 = "CD169", deg_res = summary_df_degs, genes = DEGs)[[2]]),
  annotation_name_side = "top",   
  annotation_name_rot  = 45 , annotation_name_offset = unit(2, "mm")
 )

ht_DEGs <- ComplexHeatmap::Heatmap(t(scale(t(norm_data[DEGs,]))), column_split = split.vector,
                              top_annotation = colAnn, right_annotation = rowAnn,name = "DEGs Exp.",
                              column_title = NULL, cluster_columns = F,cluster_rows = T,
                              col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
                              show_column_names = F, show_row_names = T,
                              rect_gp = grid::gpar(col = "black", lwd = 0.5), 
                              row_names_gp = grid::gpar(fontsize = 8))

lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3),
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = c("*","**","***"), type = "points", labels = c("< 0.05", "< 0.01", "< 0.001"))

pdf(paste0(output_path, "figures/DEGS_heatmap_DPvsRest.pdf"), h=10, w=7)
ComplexHeatmap::draw(ht_DEGs, column_title="DEGs", annotation_legend_list = list( lgd_pvalue,lgd_sig),
                     padding = unit(c(12, 2, 12, 2), "mm"))
dev.off()



# Comparing DP vs the other instead of DN vs the other:

rowAnn = rowAnnotation(
  `LY6E+,CD169-` = anno_simple(-log10(get_pval(group1 = "LY6E", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[1]]), 
                         col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval(group1 = "LY6E", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[2]]),
  `LY6E-,CD169+` = anno_simple(-log10(get_pval(group1 = "CD169", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[1]]), 
                           col = pvalue_col_fun,gp = gpar(col = "black"), pch = get_pval(group1 = "CD169", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[2]]),
  `LY6E+,CD169+` = anno_simple(-log10(get_pval(group1 = "DP", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[1]]), 
                            col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval(group1 = "DP", group2 = "DN", deg_res = summary_df_degs, genes = DEGs)[[2]]),
  annotation_name_side = "top",   
  annotation_name_rot  = 45 , annotation_name_offset = unit(2, "mm")
)

ht_DEGs <- ComplexHeatmap::Heatmap(t(scale(t(norm_data[DEGs,]))), column_split = split.vector,
                                   top_annotation = colAnn, right_annotation = rowAnn,name = "DEGs Exp.",
                                   column_title = NULL, cluster_columns = F,cluster_rows = T,
                                   col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
                                   show_column_names = F, show_row_names = T, 
                                   rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = 8))

lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3),
                    labels = c("1", "0.1", "0.01", "0.001"))
# and one for the significant p-values
lgd_sig = Legend(pch = c("*","**","***"), type = "points", labels = c("< 0.05", "< 0.01", "< 0.001"))

pdf(paste0(output_path, "figures/DEGS_heatmap_DNvsRest.pdf"), h=10, w=7)
# ComplexHeatmap::draw(ht, column_title="DEGs", annotation_legend_list = list( lgd_sig))
ComplexHeatmap::draw(ht_DEGs, column_title="DEGs", annotation_legend_list = list( lgd_pvalue,lgd_sig),
                     padding = unit(c(12, 2, 12, 2), "mm"))

dev.off()


# Compute signatures ------------------------------------------------------
barras.sig.df <- data.frame(readxl::read_excel(path = paste0(sinature_path, "/Clinical_signatures_Barras_Datafile_S2.xls"),sheet = 5))
head(barras.sig.df[barras.sig.df$logFC_sc<0,])
barras.sig <- barras.sig.df$Gene[!is.na(barras.sig.df$Clinical_Response_Signature_Selection)]
mulder.sig <- c("IFIT1", "MX1", "HERC5", "IFI6", "ISG15", "IFIT3", "RSAD2", "GBP1", "IFIT2", "XAF1", "PARP9", "UBE2L6", "IRF7", "PARP14", "APOL6")

mulder.sig.df <- data.frame(readxl::read_excel(path = paste0(sinature_path, "/MoMac_signatures_Mulder_1-s2.0-S1074761321002934-mmc4.xlsx"),sheet = 2))
max(mulder.sig.df$p_val_adj)
min(mulder.sig.df$avg_logFC)
mulder.sig.cxcl9 <- mulder.sig.df$gene[mulder.sig.df$cluster == 6]
mulder.sig.cxcl9

mulder.sig.isg.mono <- mulder.sig.df$gene[mulder.sig.df$cluster == 4]
mulder.sig.isg.mono

mulder.sig.isg.momac <- mulder.sig.isg.mono[mulder.sig.isg.mono%in% mulder.sig.cxcl9]

public_genesets <- list(
  "Barras_Mac_Responders" = barras.sig,
  "Mulder_Ifn_Mono" = mulder.sig.isg.mono,
  "Mulder_Mac_CXCL9" =  mulder.sig.cxcl9,
  "Mulder_Ifn_MoMac" = mulder.sig.isg.momac
)

citeSeq_mono <- read.table(citeseq_mono_markers, 
                           header = T, sep = "\t")
citeSeq_mono <- citeSeq_mono[citeSeq_mono$FDR < 0.05 & citeSeq_mono$logCPM > 0,]

htan_mac <- read.table(htan_mac_markers,sep = ";",dec = ",",header = T)
htan_mac <- htan_mac[htan_mac$FDR < 0.05 & (htan_mac$logFC > 0.8 & htan_mac$logCPM>4) & htan_mac$cluster == "Macro_CXCL9",]

atlas_genesets <- list(CD14_Ifn = citeSeq_mono$gene,
                       Macro_CXCL9 = htan_mac$gene)

all_genesets_symbols <- c(public_genesets,atlas_genesets)

# Run GSVA ----------------------------------------------------------------
allScores <- GSVA::gsva(param = GSVA::ssgseaParam(as.matrix(norm_data), all_genesets_symbols, minSize = 5))
scores_plots <- list()
for(score in names(all_genesets_symbols)){
  pData <- as.data.frame(t(allScores))
  pData$group <- factor(metadata$group, levels = c("DN", "LY6E", "CD169", "DP"))
  pData$sample <- metadata$sample
  pData$donor <- metadata$donor
  
  scores_plots[[score]] <- ggplot(pData, aes(x = group, y = .data[[score]])) +
    geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0)) +
    geom_line(aes(group = donor, color = donor), linewidth = 0.7, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar",
                 width = 0.4, fatten = 2, color = "black", linewidth = 0.5) +
    stat_compare_means(comparisons =combn(levels(factor(pData$group)), 2, simplify = FALSE),
                       method = "wilcox.test",
                       paired = T,       
                       label = "p.format")  +
    labs(title = score, x = "Group", y = "Signature") +
    theme_minimal() +
    theme()
}


pdf(paste0(output_path, "figures/signature_levels.pdf"),w=17,h=4)
scores_plots$Barras_Mac_Responders | scores_plots$Mulder_Ifn_Mono | scores_plots$Mulder_Mac_CXCL9 | scores_plots$Mulder_Ifn_MoMac
scores_plots$CD14_Ifn | scores_plots$Macro_CXCL9
dev.off()

head(allScores)
summary(colnames(allScores) == rownames(metadata))
 
# # Heatmap plots -----------------------------------------------------------
# split.vector <- factor(conditions, levels = c("DN", "LY6E", "CD169",  "DP"))
# 
# colours <- list('Condition' = descriminative_colors)
# colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
#                                             which = 'col',
#                                             col = colours)
# gsva_subset = allScores[rownames(allScores) %in% names(c(public_genesets,atlas_genesets)),, drop = F]
# 
# 
# pdf(paste0(output_path, "figures/signature_levels_heatmap.pdf"))
# 
# ht <- ComplexHeatmap::Heatmap(t(scale(t(gsva_subset))), 
#                               name = "GSVA", column_split = split.vector,
#                               top_annotation = colAnn, # right_annotation = rowAnn,
#                               column_title = NULL, cluster_columns = F,cluster_rows = T,
#                               col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
#                               show_column_names = F, show_row_names = T, rect_gp = grid::gpar(col = "black", lwd = 0.5),
#                               row_names_gp = grid::gpar(fontsize = 8))
# ComplexHeatmap::draw(ht, column_title="Signature")
# dev.off()
# 
# ComplexHeatmap::Heatmap(gsva_subset,
#                         name = "GSVA", column_split = split.vector,
#                         top_annotation = colAnn, # right_annotation = rowAnn,
#                         column_title = NULL, cluster_columns = F,cluster_rows = T,
#                         # col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
#                         show_column_names = F, show_row_names = T, rect_gp = grid::gpar(col = "black", lwd = 0.5),
#                         row_names_gp = grid::gpar(fontsize = 8))


# differential analysis on the signatures ---------------------------------

ct <- factor(metadata$group)
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)
dupcor <- duplicateCorrelation(allScores,design,block=donor)
dupcor$consensus.correlation

fit <- lmFit(allScores, design, block=donor, correlation = dupcor$consensus.correlation)

contrasts <- makeContrasts(
  CD169_vs_DN     = CD169 - DN,
  LY6E_vs_DN     =  LY6E - DN,
  DP_vs_DN =  DP - DN,
  LY6E_vs_CD169     =   LY6E - CD169,
  DP_vs_CD169     =  DP - CD169,
  DP_vs_LY6E =  DP - LY6E,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary_df_sig <- vector()
for(c in colnames(contrasts)){
  res.all = topTable(fit2, coef=c)
  res.all$group1 <- strsplit(c,"_vs_")[[1]][1]
  res.all$group2 <- strsplit(c,"_vs_")[[1]][2]
  res.all$pathway <- rownames(res.all)
  summary_df_sig <- rbind(summary_df_sig, res.all)
}
summary_df_sig_signif <- summary_df_sig[summary_df_sig$adj.P.Val <= 0.05, ]

pathways_keep <- unique(summary_df_sig_signif$pathway)
pathways_keep2 <- unique(summary_df_sig_signif$pathway[summary_df_sig_signif$group1 == "DP" & summary_df_sig_signif$group2 == "DN"])
summary(pathways_keep %in% pathways_keep2)
pathways_keep[!pathways_keep %in% pathways_keep2]

colours <- list('Condition' = descriminative_colors)
colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
                                            which = 'col',
                                            col = colours,show_annotation_name = FALSE)
# ht_sign <- ComplexHeatmap::Heatmap(t(scale(t(allScores[pathways_keep,]))), column_split = split.vector,
#                               top_annotation = colAnn, # right_annotation = rowAnn,
#                               column_title = NULL, cluster_columns = F,cluster_rows = T,
#                               col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
#                               show_column_names = F, show_row_names = T,
#                               rect_gp = grid::gpar(col = "black", lwd = 0.5),
#                               row_names_gp = grid::gpar(fontsize = 8))
# ComplexHeatmap::draw(ht_sign, column_title="Signatures")



get_pval_sig <- function(group1, group2, deg_res, genes){
  deg_res <- deg_res[which(deg_res$group1 == group1 & deg_res$group2 == group2), ]
  rownames(deg_res) <- deg_res$pathway
  pval <- deg_res[genes, "adj.P.Val"]
  names(pval) = genes
  
  is_sig = pval < 0.05
  pch = rep("*", length(pval))
  pch[!is_sig] = NA#"ns"
  pch[pval<0.01] = "**"
  pch[pval<0.001] = "***"
  
  return(list(pval, pch))
}
pvalue_col_fun = circlize::colorRamp2(c(0, 2, 3), viridis::viridis(3)) 
rowAnn_sig = rowAnnotation(
  `LY6E+,CD169-` = anno_simple(-log10(get_pval_sig(group1 = "LY6E", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[1]]), 
                         col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval_sig(group1 = "LY6E", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[2]]),
  `LY6E-,CD169+` = anno_simple(-log10(get_pval_sig(group1 = "CD169", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[1]]), 
                           col = pvalue_col_fun,gp = gpar(col = "black"), pch = get_pval_sig(group1 = "CD169", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[2]]),
  `LY6E+,CD169+` = anno_simple(-log10(get_pval_sig(group1 = "DP", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[1]]), 
                            col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval_sig(group1 = "DP", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[2]]),
  annotation_name_side = "top",   
  annotation_name_rot  = 45 , annotation_name_offset = unit(2, "mm")
)

ht_signature <- ComplexHeatmap::Heatmap(t(scale(t(allScores))), column_split = split.vector,
                              top_annotation = colAnn, right_annotation = rowAnn_sig, name = "Sig.",
                              column_title = NULL, cluster_columns = F,cluster_rows = T,
                              col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
                              show_column_names = F, show_row_names = T,
                              rect_gp = grid::gpar(col = "black", lwd = 0.5),
                              row_names_gp = grid::gpar(fontsize = 8))

lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3),
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = c("*","**","***"), type = "points", labels = c("< 0.05", "< 0.01", "< 0.001"))

pdf(paste0(output_path, "figures/Sig_heatmap.pdf"))
ComplexHeatmap::draw(ht_signature, column_title="Signatures", annotation_legend_list = list( lgd_pvalue,lgd_sig),
                     padding = unit(c(12, 2, 12, 2), "mm"))

dev.off()

rowAnn_sig = rowAnnotation(
  `LY6E+,CD169-` = anno_simple(-log10(get_pval_sig(group1 = "LY6E", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[1]]), 
                               col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval_sig(group1 = "LY6E", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[2]]),
  `LY6E-,CD169+` = anno_simple(-log10(get_pval_sig(group1 = "CD169", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[1]]), 
                               col = pvalue_col_fun,gp = gpar(col = "black"), pch = get_pval_sig(group1 = "CD169", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[2]]),
  `LY6E+,CD169+` = anno_simple(-log10(get_pval_sig(group1 = "DP", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[1]]), 
                               col = pvalue_col_fun, gp = gpar(col = "black"),pch = get_pval_sig(group1 = "DP", group2 = "DN", deg_res = summary_df_sig, genes = rownames(allScores))[[2]]),
  annotation_name_side = "top",   
  annotation_name_rot  = 45 , annotation_name_offset = unit(2, "mm"),show_annotation_name = FALSE
)

ht_signature <- ComplexHeatmap::Heatmap(t(scale(t(allScores))), column_split = split.vector,
                                        top_annotation = colAnn, right_annotation = rowAnn_sig, name = "Sig.",
                                        column_title = NULL, cluster_columns = F,cluster_rows = T,
                                        col = circlize::colorRamp2(c(-2, -1, 0 , 1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124")),
                                        show_column_names = F, show_row_names = T,
                                        rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                        row_names_gp = grid::gpar(fontsize = 10))

pdf(paste0(output_path, "figures/combined_heatmap.pdf"),w=8,h=9)

ht_all <- ht_DEGs %v% ht_signature
draw(ht_all,  annotation_legend_list = list( lgd_pvalue,lgd_sig),
     padding = unit(c(12, 2, 12, 2), "mm"))
dev.off()

summary_df_degs$group1 <- plyr::revalue(summary_df_degs$group1, replace = c("CD169" = "LY6E-,CD169+", "LY6E" = "LY6E+,CD169-", "DP" = "LY6E+,CD169+", "DN" = "LY6E-,CD169-" ))
summary_df_degs$group2 <- plyr::revalue(summary_df_degs$group2, replace = c("CD169" = "LY6E-,CD169+", "LY6E" = "LY6E+,CD169-", "DP" = "LY6E+,CD169+", "DN" = "LY6E-,CD169-" ))
write.table(summary_df_degs, file = paste0(output_path, "DEGs_res.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

