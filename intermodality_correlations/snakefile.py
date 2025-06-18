import os
import re
import pandas as pd
import scanpy as sc
import random
SEED = 448

print(config)

pbmcCiteSamples = config["pbmcCiteSeqAtlas"]["samples"]


GAMMA = config["GAMMA"]
GAMMA_seacells = [20,30,50,75,100,200]
pando_method_mc = config["pando_method_mc"]
pando_method_sc = config["pando_method_sc"]

wdir = os.getcwd()

include: "snakemakeWorkflows/pbmcMultiome/snakefile.py"
include: "snakemakeWorkflows/hspcMultiomePersad/snakefile.py"

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
    # seurat objects with activities
    "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds",
    "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
    expand("output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells"]),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"]),
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
    expand("output/pbmcMultiome/singlecells_analysis/multimodalMarkers_mainCellTypes_ttest.rds"),
    expand("output/pbmcMultiome/SuperCellMulti/g{gamma}/multimodalMarkers_mainCellTypes_surveyweightedt.rds", gamma = GAMMA),
    expand("output/pbmcMultiome/SuperCellMulti/g{gamma}/grn_object_{pandoMC}.rds", gamma = GAMMA, pandoMC = pando_method_mc),
    expand("output/pbmcMultiome/singlecells_analysis/grn_object_{pandoSC}.rds", pandoSC = pando_method_sc),
    expand("output/pbmcMultiome/SuperCellMulti/pando_metrics_summary_p_thresh0.1.RData"),
    expand("output/pbmcMultiome/SuperCellMulti/pando_metrics_summary_p_thresh0.05.RData"),
    # semi-supervised results on PBMC multiome data
    expand("output/pbmcMultiome/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv", graining = ["20","75"]),
    # generate figures
    "figures/manuscript/Figure2_pbmcMultiome_bench_v2.html",
    "figures/manuscript/figure_3_intermodality_v2.html",
    "figures/manuscript/figure_S4_intermodality.html"


# snakemake -j 30 -kps snakefile.py --configfile config/workflow.yml --use-singularity --conda-frontend conda --cluster-config config/cluster.yml --cluster "sbatch -A {cluster.account} \
#         -p {cluster.partition} \
#         -N {cluster.N} \
#         -t {cluster.time} \
#         --job-name {cluster.name} \
#         --mem {cluster.mem} \
#         --cpus-per-task {cluster.cpus-per-task}\
#         --output {cluster.output} \
#         --error {cluster.error}" --dry-run
