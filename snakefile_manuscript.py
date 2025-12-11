import os
import re
import pandas as pd
import scanpy as sc
import numpy as np
import random
SEED = 448

print(config)
GAMMA = config["GAMMA"]
GAMMA_seacell_tests = [50,75,100,200]
pando_method_mc = config["pando_method_mc"]
pando_method_sc = config["pando_method_sc"]

CiteAtlasMethod =  config["pbmcCiteSeqAtlas"]["methods"]
pbmcCiteSamples = config["pbmcCiteSeqAtlas"]["samples"]

wdir = os.getcwd()

wildcard_constraints: kernel= '(FALSE|TRUE)'

sif_file = config["sif_file"]
sif_file_edger = config["sif_file_edger"]

include: "snakemakeWorkflows/pbmcMultiome/snakefile.py"
include: "snakemakeWorkflows/bmCiteSeq/snakefile.py"
include: "snakemakeWorkflows/hspcMultiomePersad/snakefile.py"
include: "snakemakeWorkflows/pbmcCiteSeqAtlas/snakefile.py"


rule all:
  input:
    "reports/bmCiteSeq/bm_cite_analysis.html",
    "figures/manuscript/final_figures_Rmd/Figure2_bench_BMCiteseq.html",
    "figures/manuscript/final_figures_Rmd/Figure2_pbmcMultiome_bench_v2.html",
    "figures/manuscript/final_figures_Rmd/Figure4_bench_pbmc_cite_atlas_integration_v2.html",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_pairwise.txt",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_CD14.txt",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ADT_diff/edgeR_res_pairwise.txt",
    "figures/manuscript/final_figures_Rmd/figure_3_intermodality_v2.html",
    "figures/manuscript/final_figures_Rmd/figure_S4_intermodality.html"
    # expand("output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/repro.seurat.multiome.mcMetrics.rds", gamma = GAMMA_seacell_tests, rep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),


# /work/FAC/FBM/LLB/dgfeller/scrnaseq/lherault/bin/miniconda3/envs/snakemake
# snakemake -j 30 -kps snakefile_manuscript.py --configfile config/workflow_manuscript.yml --use-singularity --singularity-args "--bind /work,/scratch" --conda-frontend conda --cluster-config config/cluster.yml --cluster "sbatch -A {cluster.account} \
#         -p {cluster.partition} \
#         -N {cluster.N} \
#         -t {cluster.time} \
#         --job-name {cluster.name} \
#         {cluster.nodelist} \
#         --mem {cluster.mem} \
#         --cpus-per-task {cluster.cpus-per-task}\
#         --output {cluster.output} \
#         --error {cluster.error}" --keep-going --dry-run
