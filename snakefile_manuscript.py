import os
import re
import pandas as pd
import scanpy as sc
import random
SEED = 448

print(config)
GAMMA = config["GAMMA"]
pando_method_mc = config["pando_method_mc"]
pando_method_sc = config["pando_method_sc"]

CiteAtlasMethod =  config["pbmcCiteSeqAtlas"]["methods"]
pbmcCiteSamples = config["pbmcCiteSeqAtlas"]["samples"]

wdir = os.getcwd()

wildcard_constraints: kernel= '(FALSE|TRUE)'

sif_file = config["sif_file"]

include: "snakemakeWorkflows/pbmcMultiome/snakefile.py"
include: "snakemakeWorkflows/hspcMultiomePersad/snakefile.py"
include: "snakemakeWorkflows/bmCiteSeq/snakefile.py"
include: "snakemakeWorkflows/pbmcCiteSeqAtlas/snakefile.py"


rule all:
  input:
    "reports/bmCiteSeq/bm_cite_analysis.html",
    "figures/manuscript/Figure2_pbmcMultiome_bench_v2.html",
    "reports/pbmcCiteSeqAtlas/Correlation_MC_metrics_pbmc_cite_atlas_analysis.html",
    "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis.html"
    # "figures/manuscript/figure_3_intermodality_v2.html",
    # "figures/manuscript/figure_S4_intermodality.html",



# snakemake -j 30 -kps snakefile_manuscript.py --configfile config/workflow_manuscript.yml --use-singularity --singularity-args "--bind /work,/scratch" --conda-frontend conda --cluster-config config/cluster.yml --cluster "sbatch -A {cluster.account} \
#         -p {cluster.partition} \
#         -N {cluster.N} \
#         -t {cluster.time} \
#         --job-name {cluster.name} \
#         --mem {cluster.mem} \
#         --cpus-per-task {cluster.cpus-per-task}\
#         --output {cluster.output} \
#         --error {cluster.error}" --keep-going --dry-run
