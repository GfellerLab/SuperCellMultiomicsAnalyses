#################################################################################################
#################### BM CITE-seq dataset ##################################################
#################################################################################################
rule install_seurat_bmcite:
  output: "input/bmCiteSeq/bmcite.rds"
  singularity: config["sif_file"]
  shell: "Rscript -e 'library(Seurat);.libPaths(c(\"./input/\", .libPaths()));\
                      SeuratData::InstallData(\"bmcite\",lib = \"./input/\");\
                      bm <- SeuratData::LoadData(ds = \"bmcite\");\
                      saveRDS(bm,\"{output}\")'"

rule singlecells_rna_pca_analysis_bm_cite:
  input: dataset=  "input/bmCiteSeq/bmcite.rds"
  output: "output/bmCiteSeq/singlecells_analysis/seurat.RNA.h5ad",
  singularity: config["sif_file"]
  benchmark: "benchmark/bmCiteSeq/singlecells_analysis/seurat.RNA.txt"
  shell: "Rscript snakemakeWorkflows/bmCiteSeq/R/RNA_pca_CiteSeqCL.R -i {input.dataset} \
          -o output/bmCiteSeq/singlecells_analysis"

rule singlecells_adt_pca_analysis_bm_cite:
  input: dataset=  "input/bmCiteSeq/bmcite.rds"
  output: "output/bmCiteSeq/singlecells_analysis/seurat.ADT.h5ad",
  singularity: config["sif_file"]
  benchmark: "benchmark/bmCiteSeq/singlecells_analysis/seurat.ADT.txt"
  shell: "Rscript snakemakeWorkflows/bmCiteSeq/R/ADT_pca_CiteSeqCL.R -i {input.dataset} \
          -o output/bmCiteSeq/singlecells_analysis"

rule preprocessing_supercell_bm_cite:
  input: dataset=  "input/bmCiteSeq/bmcite.rds"
  output: "output/bmCiteSeq/singlecells_analysis/seurat_multimodal.rds",
  singularity: config["sif_file"]
  benchmark: "benchmark/bmCiteSeq/singlecells_analysis/seurat_multimodal.txt"
  shell: "Rscript snakemakeWorkflows/bmCiteSeq/R/preprocessing_seurat_for_SuperCell_CiteSeqCL.R -i {input.dataset} \
          -o output/bmCiteSeq/singlecells_analysis "

rule singlecells_wnn_analysis_bm_cite:
  input: dataset=  "input/bmCiteSeq/bmcite.rds"
  output: "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds"
  singularity: config["sif_file"]
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript snakemakeWorkflows/bmCiteSeq/R/wnnAnalysisBM_CiteCL.R -i {input.dataset} \
          -o output/bmCiteSeq/singlecells_analysis \
          -d -y {params.python} \
          -p 1:30 \
          -q 1:18"
          
rule singlecells_mofa_analysis_bm_cite:
  input: dataset=  "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
         script= "snakemakeWorkflows/bmCiteSeq/R/mofaAnalysisBM_CiteCL.R"
  output: "output/bmCiteSeq/singlecells_analysis/adata_mofa.h5ad"
  singularity: config["sif_file_mofa"]
  params: python = "/opt/conda/envs/mofa_env/bin/python"
  shell: "Rscript snakemakeWorkflows/bmCiteSeq/R/mofaAnalysisBM_CiteCL.R -i {input.dataset} \
          -o output/bmCiteSeq/singlecells_analysis \
          -y {params.python}"

rule metacell_identification_bm_cite:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seurat_multimodal.rds"
  output: "output/bmCiteSeq/SuperCellMulti/g{gamma}/SuperCellHierarchy.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  benchmark :  "benchmark/bmCiteSeq/SuperCellMulti/g{gamma}/benchIdentification.txt"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/bmCiteSeq/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:30 -q 1:18 -e TRUE -k 30 -g {wildcards.gamma} \
          -s SuperCellHierarchy"

rule data_aggregation_bm_multi:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        hierarchy = "output/bmCiteSeq/SuperCellMulti/g{gamma}/SuperCellHierarchy.rds",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/SuperCellMulti/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/SuperCellMulti/g{gamma}/metaData.csv",
          "output/bmCiteSeq/SuperCellMulti/g{gamma}/corrTablePearson.csv"

  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -m {input.hierarchy} \
          -o output/bmCiteSeq/SuperCellMulti/g{wildcards.gamma}/ \
          -g {wildcards.gamma} \
          -y SuperCell_Multi \
          -z {input.gene_protein}\
          -t {params.python}\
          -s seurat"




rule metacell_identification_bm_RNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/bmCiteSeq/SuperCellRNA/g{gamma}/SuperCellHierarchy.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  benchmark : "benchmark/bmCiteSeq/SuperCellRNA/g{gamma}/benchIdentification.txt"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/bmCiteSeq/SuperCellRNA/g{wildcards.gamma}/ \
          -p 1:30 -e TRUE -k 30 -g {wildcards.gamma} \
          -s SuperCellHierarchy"

rule data_aggregation_bm_RNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        hierarchy = "output/bmCiteSeq/SuperCellRNA/g{gamma}/SuperCellHierarchy.rds",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/SuperCellRNA/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/SuperCellRNA/g{gamma}/metaData.csv",
          "output/bmCiteSeq/SuperCellRNA/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -m {input.hierarchy} \
          -o output/bmCiteSeq/SuperCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}\
          -y SuperCell_RNA \
          -z {input.gene_protein}\
          -t {params.python}\
          -s seurat"



rule metacell_identification_bm_ADT:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.ADT.h5ad"
  output: "output/bmCiteSeq/SuperCellADT/g{gamma}/SuperCellHierarchy.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  benchmark : "benchmark/bmCiteSeq/SuperCellADT/g{gamma}/benchIdentification.txt"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/bmCiteSeq/SuperCellADT/g{wildcards.gamma}/ \
          -q 1:18 -e TRUE -k 30 -g {wildcards.gamma} \
          -s SuperCellHierarchy"

rule data_aggregation_bm_ADT:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        hierarchy = "output/bmCiteSeq/SuperCellADT/g{gamma}/SuperCellHierarchy.rds",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/SuperCellADT/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/SuperCellADT/g{gamma}/metaData.csv",
          "output/bmCiteSeq/SuperCellADT/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -m {input.hierarchy} \
          -o output/bmCiteSeq/SuperCellADT/g{wildcards.gamma}/ \
          -g {wildcards.gamma}\
          -y SuperCell_ADT \
          -z {input.gene_protein} \
          -t {params.python}\
          -s seurat"
          

rule metacell_identification_bm_metaqRNA:
 input:
       singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.RNA.h5ad"
 output: "output/bmCiteSeq/metaqRNA/g{gamma}/metaqMemberships.csv"
 singularity: config["sif_file_metaq"]
 benchmark: "benchmark/bmCiteSeq/metaqRNA/g{gamma}/benchIdentification.txt"
 shell: "/opt/conda/envs/MetaQ/bin/python python/MetaQ_CL.py --data_path {input.singlecells} \
         --data_type 'RNA' \
         --outdir output/bmCiteSeq/metaqRNA/g{wildcards.gamma}/ \
         --gamma {wildcards.gamma}"
         
rule data_aggregation_bm_metaqRNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        memberships = "output/bmCiteSeq/metaqRNA/g{gamma}/metaqMemberships.csv",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/metaqRNA/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/metaqRNA/g{gamma}/metaData.csv",
          "output/bmCiteSeq/metaqRNA/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/bmCiteSeq/metaqRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -x MetaQ-\
          -y MetaQ_RNA \
          -z {input.gene_protein} \
          -t {params.python}\
          -s seurat"

rule metacell_identification_bm_metaqADT:
  input:
      singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.ADT.h5ad"
  output: "output/bmCiteSeq/metaqADT/g{gamma}/metaqMemberships.csv"
  singularity: config["sif_file_metaq"]
  benchmark: "benchmark/bmCiteSeq/metaqADT/g{gamma}/benchIdentification.txt"
  shell: "/opt/conda/envs/MetaQ/bin/python python/MetaQ_CL.py --data_path {input.singlecells} \
         --data_type 'ADT' \
         --outdir output/bmCiteSeq/metaqADT/g{wildcards.gamma}/ \
         --gamma {wildcards.gamma}"

rule data_aggregation_bm_metaqADT:
  input:
       singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
       memberships = "output/bmCiteSeq/metaqADT/g{gamma}/metaqMemberships.csv",
       gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/metaqADT/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/metaqADT/g{gamma}/metaData.csv",
          "output/bmCiteSeq/metaqADT/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
         -o output/bmCiteSeq/metaqADT/g{wildcards.gamma}/ \
         -g {wildcards.gamma} -x MetaQ-\
         -y MetaQ_ADT \
         -z {input.gene_protein} \
         -t {params.python}\
         -s seurat"
          
rule metacell_identification_bm_metaqMulti:
  input:
      adt = "output/bmCiteSeq/singlecells_analysis/seurat.ADT.h5ad",
      rna = "output/bmCiteSeq/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/bmCiteSeq/metaqMulti/g{gamma}/metaqMemberships.csv"
  singularity: config["sif_file_metaq"]
  benchmark: "benchmark/bmCiteSeq/metaqMulti/g{gamma}/benchIdentification.txt"
  shell: "/opt/conda/envs/MetaQ/bin/python python/MetaQ_CL.py --data_path {input.rna} {input.adt} \
         --data_type 'RNA' 'ADT' \
         --outdir output/bmCiteSeq/metaqMulti/g{wildcards.gamma}/ \
         --gamma {wildcards.gamma}"

rule data_aggregation_bm_metaqMulti:
  input:
       singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
       memberships = "output/bmCiteSeq/metaqMulti/g{gamma}/metaqMemberships.csv",
       gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/metaqMulti/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/metaqMulti/g{gamma}/metaData.csv",
          "output/bmCiteSeq/metaqMulti/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
         -o output/bmCiteSeq/metaqMulti/g{wildcards.gamma}/ \
         -g {wildcards.gamma} -x MetaQ-\
         -y MetaQ_Multi \
         -z {input.gene_protein} \
         -t {params.python}\
         -s seurat"

rule metacell_identification_bm_seacellsRNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/bmCiteSeq/seacellsRNA/g{gamma}/seacellMemberships.csv"
  singularity: config["sif_file"]
  benchmark: "benchmark/bmCiteSeq/seacellsRNA/g{gamma}/benchIdentification.txt"
  shell: "python3.9 python/SEACellsCL.py -i {input.singlecells} \
          -o output/bmCiteSeq/seacellsRNA/g{wildcards.gamma}/ \
          -d 1:30 -r pca -g {wildcards.gamma}"

rule data_aggregation_bm_seacellsRNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        memberships = "output/bmCiteSeq/seacellsRNA/g{gamma}/seacellMemberships.csv",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/seacellsRNA/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/seacellsRNA/g{gamma}/metaData.csv",
          "output/bmCiteSeq/seacellsRNA/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/bmCiteSeq/seacellsRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -x SEACell-\
          -y SEACells_RNA \
          -z {input.gene_protein} \
          -t {params.python}\
          -s seurat"



rule metacell_identification_bm_seacellsMOFA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/adata_mofa.h5ad"
  output: "output/bmCiteSeq/seacellsMOFA/g{gamma}/seacellMemberships.csv"
  singularity: config["sif_file"]
  benchmark: "benchmark/bmCiteSeq/seacellsMOFA/g{gamma}/benchIdentification.txt"
  shell: "python3.9 python/SEACellsCL.py -i {input.singlecells} \
          -o output/bmCiteSeq/seacellsMOFA/g{wildcards.gamma}/ \
          -r mofa -g {wildcards.gamma}"

rule data_aggregation_bm_seacellsMOFA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        memberships = "output/bmCiteSeq/seacellsMOFA/g{gamma}/seacellMemberships.csv",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/seacellsMOFA/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/seacellsMOFA/g{gamma}/metaData.csv",
          "output/bmCiteSeq/seacellsMOFA/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/bmCiteSeq/seacellsMOFA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -x SEACell-\
          -y SEACells_MOFA \
          -z {input.gene_protein} \
          -t {params.python}\
          -s seurat"

rule metacell_identification_bm_MetaCellRNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/bmCiteSeq/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  benchmark: "benchmark/bmCiteSeq/MetaCellRNA/g{gamma}/benchIdentification.txt"
  singularity: config["sif_file"]
  shell: "python3 python/MATK_MetaCell2CL.py -i {input.singlecells} \
          -o output/bmCiteSeq/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"

rule data_aggregation_bm_MetaCellRNA:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        memberships = "output/bmCiteSeq/MetaCellRNA/g{gamma}/MetaCellMemberships.csv",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/MetaCellRNA/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/MetaCellRNA/g{gamma}/metaData.csv",
          "output/bmCiteSeq/MetaCellRNA/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/bmCiteSeq/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} \
          -y MetaCell_RNA \
          -z {input.gene_protein} \
          -t {params.python}\
          -s seurat"

rule metacell_identification_bm_seacellsADT:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seurat.ADT.h5ad",
  output: "output/bmCiteSeq/seacellsADT/g{gamma}/seacellMemberships.csv"
  singularity: config["sif_file"]
  benchmark:"benchmark/bmCiteSeq/seacellsADT/g{gamma}/benchIdentification.txt"
  shell: "python3.9 python/SEACellsCL.py -i {input.singlecells} \
          -o output/bmCiteSeq/seacellsADT/g{wildcards.gamma}/ \
          -d 1:18 -r apca -g {wildcards.gamma}"

rule data_aggregation_bm_seacellsADT:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        memberships = "output/bmCiteSeq/seacellsADT/g{gamma}/seacellMemberships.csv",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/seacellsADT/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/seacellsADT/g{gamma}/metaData.csv",
          "output/bmCiteSeq/seacellsADT/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/bmCiteSeq/seacellsADT/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -x SEACell- \
          -y SEACells_ADT \
          -z {input.gene_protein} \
          -t {params.python} \
          -s seurat"

rule random_metacell_bm_multiome:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/bmCiteSeq/gene_protein.csv"
  output: "output/bmCiteSeq/randomMetacells/g{gamma}/seurat.cite.mc.rds",
          "output/bmCiteSeq/randomMetacells/g{gamma}/metaData.csv",
           "output/bmCiteSeq/randomMetacells/g{gamma}/corrTablePearson.csv"
  singularity: config["sif_file"]
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -d TRUE \
          -o output/bmCiteSeq/randomMetacells/g{wildcards.gamma}/ \
          -g {wildcards.gamma} \
          -y random_NA \
          -z {input.gene_protein} \
          -t {params.python}\
          -s seurat"


####################################################################################################################################
################################################ Test semi-supervised  #############################################################
####################################################################################################################################



rule Semi_sup_test_bm_cite:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds"
  output: "output/bmCiteSeq/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv"
  singularity: config["sif_file"]
  params: workdir = wdir
  # benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{graining}/seurat.multiome.mc.txt"
  shell: "Rscript R/SCimplify_test_semi_sup_CL.R -i {input.singlecells} \
          -o output/bmCiteSeq/SuperCellMulti/testSemiSup/g{wildcards.graining}/ \
          -p 1:30 -q 1:18 -r RNA -a ADT -v pca -w apca -g {wildcards.graining} -s 2025 -l celltype.l2"

rule Semi_sup_test_bm_cite_rna:
  input:
        singlecells = "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds"
  output: "output/bmCiteSeq/SuperCellRNA/testSemiSup/g{graining}/results_test_semisup.csv"
  singularity: config["sif_file"]
  params: workdir = wdir
  # benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{graining}/seurat.multiome.mc.txt"
  shell: "Rscript R/SCimplify_test_semi_sup_CL.R -i {input.singlecells} \
          -o output/bmCiteSeq/SuperCellRNA/testSemiSup/g{wildcards.graining}/ \
          -p 1:30 -r RNA -v pca -g {wildcards.graining} -s 2025 -l celltype.l2"

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


rule report_bm_cite_atlas_samples:
  input:  expand("output/bmCiteSeq/{method}/g{gamma}/metaData.csv", gamma = ["20","50","75","100","200"],method = CiteAtlasMethod),
          expand("output/bmCiteSeq/{method}/g{gamma}/corrTablePearson.csv", gamma = ["20","50","75","100","200"],method = CiteAtlasMethod),
          expand("output/bmCiteSeq/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv", graining = ["20","75"]),
          expand("output/bmCiteSeq/SuperCellRNA/testSemiSup/g{graining}/results_test_semisup.csv", graining = ["20","75"])
  output: "reports/bmCiteSeq/bm_cite_analysis.html"
  shell: "touch {output}"

rule generate_figures_from_figure2_bm:
  input:
    "figures/manuscript/final_figures_Rmd/Figure2_bench_BMCiteseq.Rmd",
    expand("output/bmCiteSeq/{inputMetacells}/g{gamma}/seurat.cite.mc.rds", gamma = ["20","50","75","100","200"], inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellADT","SuperCellRNA","metaqRNA","metaqADT","metaqMulti","seacellsRNA","seacellsADT","seacellsMOFA","MetaCellRNA"]),
    expand("output/bmCiteSeq/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv", graining = ["20","75"]),
  output: "figures/manuscript/final_figures_Rmd/Figure2_bench_BMCiteseq.html"
  singularity: config["sif_file"]
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
