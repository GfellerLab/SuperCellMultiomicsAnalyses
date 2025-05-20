library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(S4Vectors)
library(patchwork)
library(getopt)


options(future.globals.maxSize= 8000*1024^2)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',     'i',1, "character", 'input RDS',
  'outdir',     'o',1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

print(opt)

combined.metacells <- readRDS(opt$input)


DepthCor(combined.metacells,reduction = "lsi")


DepthCor(combined.metacells,reduction = "integrated_lsi")



combined.metacells <- FindClusters(combined.metacells, resolution = c(c(1:7)*0.01),graph.name = "wsnn",algorithm = 3) 
DimPlot(combined.metacells,reduction = "wnn.umap", label = T,group.by = "wsnn_res.0.03")

#DimPlot(combined.metacells,reduction = "wnn.umap", label = T)


combined.metacells$log_compactness_lsi <- log(combined.metacells$compactness_lsi)
combined.metacells$log_compactness_pca <- log(combined.metacells$compactness_pca)
combined.metacells$log_INV <- log(combined.metacells$innerNormVar)


# VlnPlot(combined.metacells,features = "log_compactness_lsi")
# VlnPlot(combined.metacells,features = "log_compactness_pca")
# VlnPlot(combined.metacells,features = "log_INV")
# VlnPlot(combined.metacells,features = "size",log = T)












Idents(combined.metacells) <- "wsnn_res.0.07" 
DimPlot(combined.metacells,reduction = "wnn.umap", label = T)

# VlnPlot(combined.metacells,'nFeature_RNA')
# VlnPlot(combined.metacells,'nFeature_ATAC')



#Luminal

DefaultAssay(combined.metacells) <- "RNA"

VlnPlot(combined.metacells,features = c("Krt19","Nupr1","Lmo7","Ceacam1"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <-  RenameIdents(combined.metacells, '0' = 'Luminal_1','11' = "Luminal_2")




VlnPlot(combined.metacells,features = c("Fcgbp","Slc12a2","Lmo7","Ceacam1"),ncol = 2,pt.size = 0,assay = "RNA")



#Basal


VlnPlot(combined.metacells,features = c("Krt14","Krt15","Krt5","Krt17"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '7' = 'Basal')




# Neuroendocrine

VlnPlot(combined.metacells,features = c("Chga","Nrxn1","Hcn1","Fgf14"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '1' = 'Neuroendocrine')



# Neuroendocrine

VlnPlot(combined.metacells,features = c("Cadm2","Kcnb2","Kcnip4","Cntn4"),ncol = 2,pt.size = 0,assay = "RNA")

# Neuron

VlnPlot(combined.metacells,features = c("Plp1"),ncol = 2,pt.size = 0,assay = "RNA") + NoLegend()
FeaturePlot(combined.metacells,c("Plp1"),reduction = "wnn.umap")
combined.metacells <- RenameIdents(combined.metacells, '12' = 'Neuron')



# Basal

VlnPlot(combined.metacells,features = c("Krt14","Krt15","Krt5","Krt17"),ncol = 2,pt.size = 0,assay = "RNA")


# seminal vesicle

VlnPlot(combined.metacells,features = c("Svs5","Svs2","Cdo1","Fyb2"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '3' = 'Seminal_vesicle')



# Mesenchymal 1


VlnPlot(combined.metacells,features = c("Col1a2","Apod","Dcn","Serping1"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '2' = "Mesenchymal_1")


# Mesenchymal 2

VlnPlot(combined.metacells,features = c("Col1a2","Apod","Flt1","Serping1"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '6' = "Mesenchymal_2")


# Endothelial 1


VlnPlot(combined.metacells,features = c("Pecam1","Emcn","Selp","Flt1"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '8' = "Endothelial")



# Neutrophil


VlnPlot(combined.metacells,features = c("S100a9","Clec4d","Cxcl2","Acod1"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '4' = "Neutrophil")



# Macrophage

VlnPlot(combined.metacells,features = c("Cd74","Cd86","Cybb","Col1a"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '5' = "Macrophage")


# T cell


VlnPlot(combined.metacells,features = c("Skap1","Itk","Ptpn22","Cd3e"),ncol = 2,pt.size = 0,assay = "RNA")
combined.metacells <- RenameIdents(combined.metacells, '9' = "T_cell")


# B cell

VlnPlot(combined.metacells,features = c("Igkc","Ebf1","Mef2c","Bank1"),ncol = 2,pt.size = 0,assay = "RNA")
FeaturePlot(combined.metacells,features = c("Igkc","Ebf1","Mef2c","Bank1"))
combined.metacells <- RenameIdents(combined.metacells, '10' = "B_cell")






VlnPlot(combined.metacells,features = c("Ptprc","Epcam","Krt8"),ncol = 2,pt.size = 0,assay = "RNA")




DimPlot(combined.metacells,reduction = "wnn.umap", label = T)
DimPlot(combined.metacells,reduction = "rna.umap", label = T)
DimPlot(combined.metacells,reduction = "atac.umap", label = T)




# marker11 <- FindMarkers(combined.metacells,ident.1 = 11,only.pos = T)
# 
# marker11 <- marker11[marker11$p_val_adj< 0.05&marker11$avg_log2FC> 1,]
# marker11





my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

colors <- c('Luminal_1'='#E5D2DD',
            "Luminal_2" = '#CCC9E6',
           'Neuroendocrine'='#53A85F',
           'Basal' = '#F1BB72',
           "Seminal_vesicle"= '#F3B1A0',
           'Mesenchymal_1'='#D6E7A3',
           'Mesenchymal_2'='#57C3F3',
           "Endothelial"='#476D87',
           "Neutrophil"='#E95C59',
           "Macrophage"='#AB3282',
           "T_cell" = '#23452F',
           "B_cell" = '#BD956A',
           "Neuron" = '#8C549C')

#####
# Idents(combined.metacells) <- combined.metacells$wsnn_res.0.1
# 
#  combined.metacells <- RenameIdents(combined.metacells, '0' = 'Luminal_1','1' = 'Mesenchymal_1','2' = 'Neuroendocrine_1','3' = 'Seminal_vesicle','4' = 'Neutrophil','5' = 'Macrophage', '6' = 'Mesenchymal_2','7' = 'Basal','8' = 'Endothelial','9' = 'T_cell','10' = 'Neuroendocrine_2','11' = 'B_cell','12' = 'Neuron')

combined.metacells@meta.data$celltype <- factor(Idents(combined.metacells),                                                levels = c("Basal","Luminal_1","Luminal_2","Neuroendocrine",
                                                                                                                                                 "Seminal_vesicle","Mesenchymal_1","Mesenchymal_2",
                                                                                                                                                 "Endothelial","Neutrophil","Macrophage",
                                                                                                                                                 "T_cell","B_cell","Neuron"))

Idents(combined.metacells) <- "celltype"


pdf(paste0(opt$outdir,"/allAnnotatedMetacells_wnn.umap.pdf"))
plot1 <- DimPlot(combined.metacells,reduction = "wnn.umap",cols = colors,label = T,repel = T)
plot1 + NoLegend()
#LabelClusters(plot1, id = "ident", color = unique(ggplot_build(plot1)$data[[1]]$colour), size = 5, repel = T,  box.padding = 3)
dev.off()

ggplot(combined.metacells@meta.data,aes(x= celltype,y=log(compactness_pca),fill = celltype)) + geom_boxplot() +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggplot(combined.metacells@meta.data,aes(x= celltype,y=log(compactness_lsi),fill = celltype)) + geom_boxplot() +
  scale_x_discrete(guide = guide_axis(angle = 45))






combined.metacells@meta.data$sample <- factor(combined.metacells$orig.ident,
                                                  levels = c("WT_1","WT_2",
                                                             "2W_1","2W_2",
                                                             "1M_1","1M_2",
                                                             "2_5M_1","2_5M_2","2_5M_3","2_5M_4",
                                                             "3_5M_1","3_5M_2",
                                                             "4_5M_1","4_5M_2",
                                                             "6M"))

combined.metacells$time <- combined.metacells$orig.ident
combined.metacells$time[startsWith(combined.metacells$orig.ident,"1M")] <- "1M"
combined.metacells$time[startsWith(combined.metacells$orig.ident,"2_5M")] <- "2_5M"
combined.metacells$time[startsWith(combined.metacells$orig.ident,"2W")] <- "2W"
combined.metacells$time[startsWith(combined.metacells$orig.ident,"3_5M")] <- "3_5M"
combined.metacells$time[startsWith(combined.metacells$orig.ident,"4_5M")] <- "4_5M"
combined.metacells$time[startsWith(combined.metacells$orig.ident,"6M")] <- "6M"
combined.metacells$time[startsWith(combined.metacells$orig.ident,"WT")] <- "WT"

combined.metacells@meta.data$time <- factor(combined.metacells$time,levels = c("WT","2W","1M",
                                                                     "2_5M","3_5M","4_5M",
                                                                     "6M"))



ggplot(combined.metacells@meta.data,aes(x=sample,fill = time)) + geom_bar(stat = "count")  +
  scale_x_discrete(guide = guide_axis(angle = 45)) + ylab("metacell counts")







immune_types <- c("T_cell","B_cell","Neutrophil","Macrophage")

ggplot(combined.metacells@meta.data[combined.metacells$celltype %in% immune_types,],aes(x=sample,fill = time)) + geom_bar(stat = "count")  +
  scale_x_discrete(guide = guide_axis(angle = 45)) + ylab("immune metacell counts")


combined.metacells$major_type <- combined.metacells$celltype
combined.metacells$major_type <- as.character(gsub(x=combined.metacells$major_type, pattern = "_[0-9]",replacement = ""))
# combined.metacells$major_type[combined.metacells$major_type == "Mesenchymal"] <- as.character(combined.metacells$major_type_res.0.2[combined.metacells$major_type == "Mesenchymal"])

DimPlot(combined.metacells,reduction = "wnn.umap",label = T,repel = T,group.by = "major_type")

colors <- c('Luminal'='#E5D2DD',
           'Neuroendocrine'='#53A85F',
           'Basal' = '#F1BB72',
           "Seminal_vesicle"= '#F3B1A0',
           'Mesenchymal'='#D6E7A3',
           # 'Mesenchymal_2'='#57C3F3',
           # 'Mesenchymal_3'='#91D0BE',
           "Endothelial"='#476D87',
           "Neutrophil"='#E95C59',
           "Macrophage"='#AB3282',
           "T_cell" = '#23452F',
           "B_cell" = '#BD956A',
           "Neuron" = '#8C549C')

combined.metacells@meta.data$major_type <- factor(combined.metacells$major_type,                                                
                                                  levels = (c("Luminal","Neuroendocrine","Basal",
                                                             "Seminal_vesicle","Mesenchymal","Endothelial",
                                                             "Neutrophil","Macrophage",
                                                             "T_cell","B_cell","Neuron"))) 

pdf(paste0(opt$outdir,"/allMetacells_majorTypes_wnn.umap.pdf"))
DimPlot(combined.metacells,reduction = "wnn.umap",group.by = "major_type",cols = colors)

dev.off()


## cell type abundances analyses
library(reshape2)
#combined.metacells$major_type <- droplevels(combined.metacells$major_type)
smpCounts <- aggregate(combined.metacells$size, by=list(sample = combined.metacells$sample,
                                                        major_type = combined.metacells$major_type,
                                                        celltype = combined.metacells$celltype,
                                                        time = combined.metacells$time), FUN=sum)
ggplot(smpCounts,aes(x=sample,y = x,fill = time)) + geom_bar(stat = "identity")  +
  scale_x_discrete(guide = guide_axis(angle = 45)) + ylab("sc counts")





library(colorspace)
darker_colors <- darken(colors, 0.2)

contingencyTable <- xtabs(x ~ major_type+time,data = smpCounts)
colSums(contingencyTable)

freqMatrix <- apply(contingencyTable,1,FUN = function(x){x/colSums(contingencyTable)})
# res <- chisq.test(contingencyTable)
# Roe <- res$observed/res$expected

freqMatrix_df <- melt(freqMatrix)

freqMatrix_df$major_type <- factor(freqMatrix_df$major_type ,
                                   levels = c("Basal","Endothelial","Luminal","Neuroendocrine","B_cell","T_cell","Macrophage","Neutrophil","Seminal_vesicle","Mesenchymal","Neuron"))

freqMatrix_df$group <- freqMatrix_df$major_type

freqMatrix_df$group[freqMatrix_df$time == "WT"] <- NA

pdf(paste0(opt$outdir,"/celltype_abundances_plots.pdf"))
ggplot(freqMatrix_df,aes(x = time,y=value*100,color=major_type,group = group)) + geom_point() + geom_line() + facet_wrap("~major_type") +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + scale_color_manual(values = darker_colors) + ylab("% total cells (mean of samples)")



ggplot(freqMatrix_df[freqMatrix_df$major_type %in% c("Luminal","Neuroendocrine"),],aes(x = time,y=value*100,color=major_type,group = group)) + geom_point() + geom_line() +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + scale_color_manual(values = darker_colors) + ylab("% total cells")



ggplot(freqMatrix_df[freqMatrix_df$major_type %in% c("Luminal","Neuroendocrine","Neutrophil"),],aes(x = time,y=value*100,color=major_type,group = group)) + geom_point() + geom_line() +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + scale_color_manual(values = darker_colors) + ylab("% total cells")

dev.off()




DefaultAssay(combined.metacells) <- "RNA"
genes <- c("Krt8","Chga","Krt5","Svs5","Col1a2","Cdh5","S100a9","C1qa","Cd3e","Cd19","Plp1")
genes

combined.metacells@meta.data$major_type <- factor(combined.metacells$major_type,                                                
                                                  levels = rev(c("Luminal","Neuroendocrine","Basal",
                                                             "Seminal_vesicle","Mesenchymal","Endothelial",
                                                             "Neutrophil","Macrophage",
                                                             "T_cell","B_cell","Neuron"))) 

pdf(paste0(opt$outdir,'/vln_stack_all_clusters.pdf'))
VlnPlot(object = combined.metacells,features = genes,stack = T,group.by = "major_type",fill.by= "ident",cols = colors,same.y.lims = F)  + NoLegend() + theme(axis.title.y = element_blank())

dev.off()


# 
# DefaultAssay(combined.metacells) <- "ATAC"
# frags <- Fragments(combined.metacells)  # get list of fragment objects
# Fragments(combined.metacells) <- NULL  # remove fragment information from assay
# 
# new.paths <- list.files('../../output/10xMultiomeProstateHan22/logNorm/',
#                         pattern='MC_atac_fragments.tsv.gz',
#                         recursive = T,
#                         full.names = T)
# 
# new.paths <- new.paths[grepl(x = new.paths,pattern = "aggregated_fragment_file/MC_atac_fragments.tsv.gz")]
# new.paths <- new.paths[!grepl(x = new.paths,pattern = "tbi")]
# 
# 
# for (i in seq_along(frags)) {
#   print(i)
#   frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[i]) # update path
# }
# 
# #Error!!!
# Fragments(combined.metacells) <- frags # assign updated list back to the object

write.csv(combined.metacells@meta.data,paste0(opt$outdir,"/combined.metacells_all_anno.csv"))
write.csv(combined.metacells@reductions$wnn.umap@cell.embeddings,paste0(opt$outdir,"/combined.metacells_wnn.umap_metacells_embed.csv"))

## Coverage plots
DefaultAssay(combined.metacells) <- "ATAC"
genes <- c("Krt8","Chga","Krt5","Svs5","Col1a2","Cdh5","S100a9","C1qa","Cd3e","Cd19","Plp1")
combined.metacells@meta.data$major_type <- factor(combined.metacells$major_type,                                                
                                                  levels = c("Luminal","Neuroendocrine","Basal",
                                                             "Seminal_vesicle","Mesenchymal","Endothelial",
                                                             "Neutrophil","Macrophage",
                                                             "T_cell","B_cell","Neuron"))

#genes <- c("S100a9","C1qa","Cd3e","Cd19")

TSS <- GetTSSPositions(Annotation(combined.metacells))
TSS <- TSS[TSS$gene_name%in% genes,]

regions <- lapply(X = genes, function(X){
  TSSgene <- TSS[TSS$gene_name == X]
  region <- Extend(TSSgene,upstream = 100,downstream = 100)
  return(region)
  })

regions

p <- CoveragePlot(
    object = combined.metacells,
    region = regions,ncol = length(TSS),
    peaks = FALSE,color = colors,
  ) & scale_fill_manual(values = colors) 

p$patches$plots[[1]] <- p$patches$plots[[1]]  & theme(
                             axis.text.x = element_blank(),
                             #axis.title.x = element_blank(),
                             axis.ticks.x= element_blank(),axis.ticks.y = element_blank()) & xlab(gsub(x=p$patches$plots[[1]]$labels$x,pattern = " position \\(bp\\)",replacement = ""))

for (i in c(2:length(TSS))) {
  p$patches$plots[[i]] <- p$patches$plots[[i]] & theme(strip.text.y.left = element_blank(),strip.background = element_blank())  & theme(strip.text.y.left = element_blank(),
                             strip.background = element_blank(),
                             axis.title.y = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x= element_blank(),axis.ticks.y = element_blank())& xlab(gsub(x=p$patches$plots[[i]]$labels$x,pattern = " position \\(bp\\)",replacement = ""))
 

}

pdf(paste0(opt$outdir,"/coverage_plots_all_cells.pdf"))
p
dev.off()


#save immune

combined.metacells[,combined.metacells$major_type %in% c("T_cell",'B_cell',"Macrophage","Neutrophil")]











Idents(combined.metacells) <- combined.metacells@meta.data$major_type 
  
genes <- c("Krt8","Chga","Krt5","Svs5","Col1a2","Cdh5","S100a9","C1qa","Cd3e","Cd19","Plp1")

genes <- c("S100a9","C1qa","Cd3e","Cd19")

TSS <- GetTSSPositions(Annotation(combined.metacells))
TSS <- TSS[TSS$gene_name%in% genes,]

regions <- lapply(X = genes, function(X){
  TSSgene <- TSS[TSS$gene_name == X]
  region <- Extend(TSSgene,upstream = 500,downstream = 100)
  return(region)
  })

regions

p <- CoveragePlot(
    object = combined.metacells,
    region = regions,ncol = length(TSS),
    peaks = F,color = colors, 
  ) & scale_fill_manual(values = colors) 

for (i in c(2:length(TSS))) {
  p$patches$plots[[i]] <- p$patches$plots[[i]] & theme(strip.text.y.left = element_blank(),strip.background = element_blank())

}


p$patches$plots[[i]] & theme(strip.text.y.left = element_blank(),
                             strip.background = element_blank(),
                             title = element_blank(),
                             axis.text.x = element_blank())

pdf(paste0(opt$outdir,"/coverage_plot_immune_genes_all_cells.pdf"))

p

dev.off()










# genes <- c("Krt8","Chga","Krt5","Svs5","Col1a2","Cdh5","S100a9","C1qa","Cd3e","Cd19","Plp1")
# DefaultAssay(combined.metacells) <- "ATAC"
# 
png(paste0(opt$outdir,'CoveragePlots.png'),height = 500,width = 1200)
 plot(CoveragePlot(
    object = combined.metacells,
    region = c("S100a9","C1qa","Cd3e","Cd19"),
    ncol = 4,
    # extend.upstream = -3000,
    # extend.downstream = 1000,
    peaks = FALSE
  ))+ scale_fill_manual(values = colors)
dev.off()


DefaultAssay(combined.metacells) <- "ATAC"


for (g in  genes) {
  print(g)
  png(paste0(opt$outdir,'CoveragePlots_',g,'.png'), width = 380, height = 280, units = "px")
  plot <- CoveragePlot(
    object = combined.metacells,
    region = g,
    expression.assay = "RNA",
    # extend.upstream = -8000,
    # extend.downstream = 1000,
    peaks = FALSE
  )+ scale_fill_manual(values = colors)
  plot(plot)
  dev.off()
}

saveRDS(combined.metacells[,combined.metacells$major_type %in% c("T_cell",'B_cell',"Macrophage","Neutrophil")],
        paste0(opt$outdir,"/immune.combined.metacells.rds"))





