import os
import re
import pandas as pd
import scanpy as sc
import random
SEED = 448

print(config)

#baseDir = config["dir"]["base"]
#inputDataDir = config["dir"]["inputData"]

# cellTypesMixology = config["testKernel"]["cellTypesMixology"]
#
# cellTypesBmCiteSeq = config["testKernel"]["cellTypesBmCiteSeq"]
# cellTypesPBMCcellType10Xmultiome = config["testKernel"]["cellTypesPBMCcellType10Xmultiome"]

# seuratDatasets = ["bmcite","pbmcMultiome"]

# GAMMA_CITE_ATLAS = config["pbmcCiteSeqAtlas"]["GAMMA"]
#
pbmcCiteSamples = config["pbmcCiteSeqAtlas"]["samples"]

# CiteAtlasMethod =  config["pbmcCiteSeqAtlas"]["methods"]
#
# velocitySamples = config["multiomicVelocity"]["samples"]
#
# multiveloRuns = config["multiomicVelocity"]["runs"]

GAMMA = config["GAMMA"]
GAMMA_seacells = [20,30,50,75,100,200]
pando_method_mc = config["pando_method_mc"]
pando_method_sc = config["pando_method_sc"]

wdir = os.getcwd()

#wildcard_constraints: pbmcCiteSmp= '^P[1-9]'
# wildcard_constraints: kernel= '(FALSE|TRUE)'


include: "snakemakeWorkflows/pbmcMultiome/snakefile.py"
include: "snakemakeWorkflows/hspcMultiomePersad/snakefile.py"
include: "snakemakeWorkflows/pbmcCiteSeqAtlas/snakefile.py"

rule all:
  input:
    "input/pbmcMultiome/pbmcMultiome_broad_peaks.rds",
    "output/hspcMultiomePersad/singlecells_analysis/seurat.RNA.h5ad",
    "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad",
    "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
    "output/hspcMultiomePersad/singlecells_analysis/seurat_multimodal.rds",
    "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
    "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds",
    "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
    # metacell outputs
    expand("output/hspcMultiomePersad/SuperCellMulti/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    # expand("output/pbmcMultiome/SuperCellATAC/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    # expand("output/pbmcMultiome/MetaCellRNA/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    # expand("output/pbmcMultiome/seacellsRNA/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    # expand("output/pbmcMultiome/seacellsATAC/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    # seurat objects with activities
    "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds",
    "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
    # expand("output/hspcMultiomePersad/SuperCellMulti/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA),
    expand("output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells"]),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    # #random metacellls output
    # expand("output/hspcMultiomePersad/randomMetacells/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    # expand("output/pbmcMultiome/randomMetacells/g{gamma}/seurat.multiome.mc.rds", gamma = GAMMA),
    # metacell metrics
    expand("output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.mcMetrics.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells"]),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mcMetrics.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    # multimodal markers
    expand("output/pbmcMultiome/singlecells_analysis/multimodalMarkers_ttest.rds"),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/multimodalMarkers_surveyweightedt.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    expand("output/pbmcMultiome/singlecells_analysis/CorrTables.rds"),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/CorrTables.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    expand("output/hspcMultiomePersad/singlecells_analysis/CorrTables.rds"),
    expand("output/hspcMultiomePersad/{inputMetacells}/g{gamma}/CorrTables.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells"]),
    # expand("output/hspcMultiomePersad/randomMetacells/g{gamma}/CorrTables.rds", gamma = GAMMA),
    expand("output/pbmcMultiome/singlecells_analysis/multimodalMarkers_mainCellTypes_ttest.rds"),
    expand("output/pbmcMultiome/SuperCellMulti/g{gamma}/multimodalMarkers_mainCellTypes_surveyweightedt.rds", gamma = GAMMA),
    # expand("output/hspcMultiomePersad/singlecells_analysis/multimodalMarkers_ttest.rds"),
    # expand("output/hspcMultiomePersad/SuperCellMulti/g{gamma}/multimodalMarkers_surveyweightedt.rds", gamma = GAMMA),
    #
    expand("output/pbmcMultiome/SuperCellMulti/g{gamma}/grn_object_{pandoMC}.rds", gamma = GAMMA, pandoMC = pando_method_mc),
    expand("output/pbmcMultiome/singlecells_analysis/grn_object_{pandoSC}.rds", pandoSC = pando_method_sc),
    expand("output/pbmcMultiome/SuperCellMulti/pando_metrics_summary_p_thresh0.1.RData"),
    expand("output/pbmcMultiome/SuperCellMulti/pando_metrics_summary_p_thresh0.05.RData")
    # "reports/pbmcMultiome/pbmcMultiome.html"
    # "input/GSM5008737_ADT_3P/features.tsv.gz",
    # expand("output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat_multimodal.rds", pbmcCiteSmp = pbmcCiteSamples)
    # "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis.html",
    #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/seuratCombinedWNN_2.rds",
    #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scATOMIC/seuratCombinedWNN.rds",
    #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scGate/seuratCombinedWNN.rds",
    #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_unsup/seuratCombinedWNN.rds"
    #"reports/pbmcMultiome/pbmcMultiome.html"
    #"reports/hspcMultiomePersad/hspcMultiomePersad.html"
    #"reports/pbmcCiteSeqAtlas/Correlation_MC_metrics_pbmc_cite_atlas_analysis.html"

# rule all:
#   input: expand("output/multiomicVelocity/{velocitySmp}/{runMv}/multivelo_result.h5ad",velocitySmp = velocitySamples,runMv = multiveloRuns),
#         #"reports/testKernel/kernel_analyses.html",
#          expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTable.csv",gamma = GAMMA,inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
#          expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTablePearson.csv",gamma = GAMMA,inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
#          expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/chromVarCorrTable.csv",gamma = GAMMA,inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
#          #"output/covidRen21/supervisedMajorTypeMetacellsStacasIntegration/all_cells_integrated.stacas.Rds",
#          "output/covidRen21/supervisedMajorTypeMetacellsStacasIntegration/myelo_cells_integrated.stacas.Rds",
#          "output/pbmcMultiome/singlecells_analysis/corrTable.csv",
#          "output/pbmcMultiome/singlecells_analysis/corrTablePearson.csv",
#          "output/pbmcMultiome/singlecells_analysis/chromVarCorrTable.csv",
#          "output/pbmcMultiome/singlecells_analysis/filtered.exp.scenic.csv",
#          "output/PBMC_Boukhaled/metacells/call_peaks_on_integrated_rna/mapped.multiomics.metacells.macs2.peaks.rds",
#          "output/pbmcMultiome/SuperCellMulti/g10/seurat.multiome.mc.rds",
#          #"output/pbmcMultiome/SuperCellMulti/g10/SCENIC/cis_target/regulons.csv",
#          # "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis.html",
#          "reports/pbmcCiteSeqAtlas/Correlation_MC_metrics_pbmc_cite_atlas_analysis.html",
#          "reports/bmCiteSeq/bm_cite_analysis.html",
#          "input/prostateCancer10xMultiome/rna.All.combined.integrated.Rdata",
#          "input/prostateCancer10xMultiome/atac.prostate.merge.Rdata",
#          "input/prostateCancer10xMultiome/wnn.all.combined.integrated_meta.data.csv",
#          #"output/10xMultiomeProstateHan22/logNorm/combined.metacells.rds",
#          #"output/10xMultiomeProstateHan22/logNorm/macs2/combined.metacells.rds",
#          #"output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/combined.metacells.rds",
#          #"output/10xMultiomeProstateHan22/SCT/combined.metacells.rds",
#          #"output/10xMultiomeProstateHan22/gamma_5/logNorm/combined.metacells.rds",
#          "output/covidRen21/supervisedMajorTypeMetacellsSeuratIntegration/integrated.seurat.rpca.Rds",
#          "output/covidRen21/unsupMetacellsSeuratIntegration/integrated.seurat.rpca.Rds",
#           #"output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune_hoccomoco_motifs/combined.metacells.with.macs.peak.rds",
#           #"output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune_JASPAR_motifs/combined.metacells.with.macs.peak.rds",
#          #"output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune.combined.metacells.rds",
#          #"output/pbmcMultiome/singlecells_analysis/SCENIC/AUCell/regulons_enrichment.csv",
#          #"output/pbmcMultiome/SuperCellMulti/g10/SCENIC/AUCell/regulons_enrichment.csv",
#           #"output/pbmcMultiome/SuperCellMulti/g10/SCENIC/cis_target/regulons.json",
#           #"output/pbmcMultiome/singlecells_analysis/SCENIC/cis_target/regulons.json",
#           "reports/pbmcMultiome/bench_global_trans_reg.html",
#            "reports/pbmcMultiome/bench_monocytes_trans_reg.html"
