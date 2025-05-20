import os
import re
import pandas as pd
import scanpy as sc
import random
SEED = 448

print(config)

#baseDir = config["dir"]["base"]
#inputDataDir = config["dir"]["inputData"]

SuperCellMultiomics = config["R4.1_env"]
seuratV5 = config["seuratV5_env"]
han_env = config["han_env"]
linger_env = config["linger_env"]
metacell2_0_9_env = config["metacell2_0_9_env"]
scenic_env = config["scenic_env"]
SuperCellMultiomics_pyenv = config["seaMetaCells_env"]
JASPAR2024 = config["JASPAR2024_env"]
SupervisedTools = config["scATOMIC_env"]
sif_file = config["sif_file"]

## Retrieving python env for computing compactness/velocity in R with reticulate
#Seacells = config["seacells"]
snakemake_conda_env = [f for f in os.listdir('.snakemake/conda') if f.endswith(".yaml")]
for e in snakemake_conda_env:
    with open(".snakemake/conda/"+e) as f:
        first_line = f.readline().strip('\n')
        if "seaMetaCells" in first_line:
            Seacells = ".snakemake/conda/"+os.path.splitext(e)[0]
            break
        
print(Seacells)


cellTypesMixology = config["testKernel"]["cellTypesMixology"]

cellTypesBmCiteSeq = config["testKernel"]["cellTypesBmCiteSeq"]
cellTypesPBMCcellType10Xmultiome = config["testKernel"]["cellTypesPBMCcellType10Xmultiome"]

seuratDatasets = ["bmcite","pbmcMultiome"]

GAMMA = config["GAMMA"]

GAMMA_CITE_ATLAS = config["pbmcCiteSeqAtlas"]["GAMMA"]

pbmcCiteSamples = config["pbmcCiteSeqAtlas"]["samples"]

CiteAtlasMethod =  config["pbmcCiteSeqAtlas"]["methods"]

velocitySamples = config["multiomicVelocity"]["samples"] 

multiveloRuns = config["multiomicVelocity"]["runs"] 

wdir = os.getcwd() 

#wildcard_constraints: pbmcCiteSmp= '^P[1-9]'
wildcard_constraints: kernel= '(FALSE|TRUE)'


include: "snakemakeWorkflows/installations/snakefile.py"
include: "snakemakeWorkflows/10xMultiomeProstateHan22/snakefile.py"
include: "snakemakeWorkflows/multiomicVelocity/snakefile.py"
include: "snakemakeWorkflows/testKernel/snakefile.py"
include: "snakemakeWorkflows/pbmcMultiome/snakefile.py"
include: "snakemakeWorkflows/hspcMultiomePersad/snakefile.py"
include: "snakemakeWorkflows/pbmcCiteSeqAtlas/snakefile.py"
include: "snakemakeWorkflows/bmCiteSeq/snakefile.py"
include: "snakemakeWorkflows/covidRen21/snakefile.py"
include: "snakemakeWorkflows/PBMC_Boukhaled/snakefile.py"



rule all:
  input:  "reports/bmCiteSeq/bm_cite_analysis.html",
          "reports/pbmcMultiome/pbmcMultiome.html",
          "reports/pbmcCiteSeqAtlas/Correlation_MC_metrics_pbmc_cite_atlas_analysis.html",
           "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis.html"
          # "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis.html",
          #  "output/testKernel/pbmcMultiome/detectRarePop/all_seeds_processed_at_all_gammas"
        #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/seuratCombinedWNN_2.rds",
        #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scATOMIC/seuratCombinedWNN.rds",
         #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scGate/seuratCombinedWNN.rds",
        #"output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_unsup/seuratCombinedWNN.rds"
    #"reports/bmCiteSeq/bm_cite_analysis.html"
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



