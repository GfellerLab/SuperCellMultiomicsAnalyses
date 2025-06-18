#################################################################################################
#################### hspc 10X multiome dataset ##################################################
#################################################################################################


urlsPersad = list(config["urlHspcMultiomePersad"].values())


rule download_data_from_zenodo:
  output: "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz",
          "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz.tbi",
          "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz",
          "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz.tbi",
          "input/hspcMultiomePersad/cd34_multiome_atac.h5ad",
          "input/hspcMultiomePersad/cd34_multiome_rna.h5ad"
  params: file_urls = urlsPersad
  shell: "cd input/hspcMultiomePersad/;\
          wget {params.file_urls}"

rule singlecells_rna_pca_analysis_hspc_multiome:
  input:
    adata = "input/hspcMultiomePersad/cd34_multiome_rna.h5ad"
  output: "output/hspcMultiomePersad/singlecells_analysis/seurat.RNA.h5ad"
  singularity: "config/supercell_multiomics.sif"
  benchmark: "benchmark/hspcMultiomePersad/singlecells_analysis/seurat.RNA.txt"
  shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/RNA_pca_hspcMultiomeCL.R -i {input.adata} \
          -o output/hspcMultiomePersad/singlecells_analysis"

rule preprocessing_supercell_hspc_multiome:
  input: fragment_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz",
         fragment_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz",
         index_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz.tbi",
         index_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz.tbi",
         atac = "input/hspcMultiomePersad/cd34_multiome_atac.h5ad",
         rna = "input/hspcMultiomePersad/cd34_multiome_rna.h5ad"
  output: "output/hspcMultiomePersad/singlecells_analysis/seurat_multimodal.rds"
  singularity: "config/supercell_multiomics.sif"
  benchmark: "benchmark/hspcMultiomePersad/singlecells_analysis/seurat_multimodal.txt"
  shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/preprocessing_seurat_for_SuperCell_hspcMultiomeCL.R -i {input.rna}\
          -j {input.atac} -f {input.fragment_rep1}+{input.fragment_rep2} \
          -o output/hspcMultiomePersad/singlecells_analysis"

rule wnn_single_cell_hspc_multiome:
  input: fragment_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz",
         fragment_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz",
         index_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz.tbi",
         index_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz.tbi",
         atac = "input/hspcMultiomePersad/cd34_multiome_atac.h5ad",
         rna = "input/hspcMultiomePersad/cd34_multiome_rna.h5ad"
  output: "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds","output/hspcMultiomePersad/singlecells_analysis/selectedGenes.txt"
  singularity: "config/supercell_multiomics.sif"
  benchmark: "benchmark/hspcMultiomePersad/singlecells_analysis/seuratWNN.txt"
  resources:
        mem_mb = 100000
  shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/wnnAnalysis10xMultiomeCL.R -i {input.rna}\
          -j {input.atac} -f {input.fragment_rep1}+{input.fragment_rep2} \
          -p 1:50 -q 2:50 \
          -o output/hspcMultiomePersad/singlecells_analysis -t {threads}"

rule metacell_identification_hspc_multiome:
  input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat_multimodal.rds"
  output: "output/hspcMultiomePersad/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  singularity: "config/supercell_multiomics.sif"
  params: workdir = wdir
  benchmark :  "benchmark/hspcMultiomePersad/SuperCellMulti/g{gamma}/seurat.multiome.mc.txt"
  resources:
        mem_mb = 100000
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} \
          -o output/hspcMultiomePersad/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:50 -q 2:50 -e TRUE -k 30 -g {wildcards.gamma} -b"

rule data_aggregation_hspc_SuperCellMulti:
  input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  output: "output/hspcMultiomePersad/SuperCellMulti/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/supercell_multiomics.sif"
  params: workdir = wdir
  resources:
        mem_mb = 100000
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/SuperCellMulti/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

rule random_metacell_hspc_multiome:
  input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds"
  output: "output/hspcMultiomePersad/randomMetacells/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/supercell_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -d TRUE \
          -o output/hspcMultiomePersad/randomMetacells/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

rule GA_computation_sc_hspc:
  input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
  output: "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds"
  singularity: "config/supercell_multiomics.sif"
  params: genome = config["hspcMultiomePersad"]["genome"]
  resources:
        mem_mb = 100000
  shell: "Rscript R/compute_activities.R -i {input.singlecells} \
          -o output/hspcMultiomePersad/singlecells_analysis/ \
          -w {threads} -g {params.genome} -u 2000 -d 100 -a ATAC -r -c -p {input.motifs}"

rule GA_computation_mc_hspc:
  input:
        metacells = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.mc.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
  output: "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds"
  singularity: "config/supercell_multiomics.sif"
  params: genome = config["hspcMultiomePersad"]["genome"]
  resources:
        mem_mb = 100000
  shell: "Rscript R/compute_activities.R -i {input.metacells} \
          -o output/hspcMultiomePersad/{wildcards.inputMetacells}/g{wildcards.gamma}/ \
          -w {threads} -g {params.genome} -u 2000 -d 100 -a ATAC -r -c -p {input.motifs}"


rule mc_metrics_computation_hspc:
  input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds",
        metacells = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds"
  output: "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.mcMetrics.rds"
  singularity: "config/supercell_multiomics.sif"
  params:
        outdir = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/",
        python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python"
  shell: "Rscript R/compute_MC_metrics.R -s {input.singlecells} -o {params.outdir} -m {input.metacells} -n 1:50 -q 2:50 -e {params.python}"


####################################################################################################################################
########################################## Benchmark and Correlations  #############################################################
####################################################################################################################################

rule atacRnaSinglecellAnalysis_hspc:
    input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
    output:
        "output/hspcMultiomePersad/singlecells_analysis/wilcox_rna.rds",
        "output/hspcMultiomePersad/singlecells_analysis/wilcox_motifs.rds",
        "output/hspcMultiomePersad/singlecells_analysis/multimodalMarkers_ttest.rds",
        "output/hspcMultiomePersad/singlecells_analysis/multimodalMarkers_wilcoxonAUC.rds",
        "output/hspcMultiomePersad/singlecells_analysis/CorrTables.rds"
    params:
        outdir = "output/hspcMultiomePersad/singlecells_analysis/",
        genome = config["hspcMultiomePersad"]["genome"],
        python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python"
    resources:
        mem_mb = 256000
    singularity: "config/supercell_multiomics.sif"
    shell: "Rscript R/AtacRnaCorrAnalysis.R -i {input.singlecells} -o {params.outdir} -s {input.singlecells}"

rule atacRnaMetacellAnalysis_hspc:
    input:
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds",
        metacells = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
    output:
        "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/multimodalMarkers_surveyweightedt.rds",
        "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/CorrTables.rds"
    params:
        outdir = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/",
        genome = config["hspcMultiomePersad"]["genome"]
    resources:
        mem_mb = 100000
    singularity: "config/supercell_multiomics.sif"
    shell: "Rscript R/AtacRnaCorrAnalysis.R -i {input.metacells} -o {params.outdir} -s {input.singlecells}"

####################################################################################################################################
########################################## Figures #################################################################################
####################################################################################################################################

rule generate_figures_from_figureS4:
  input: "figures/manuscript/figure_S4_intermodality.Rmd",
        "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.activities.rds",
        expand("output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti"])
  output: "figures/manuscript/figure_S4_intermodality.html"
  singularity: "config/supercell_multiomics_v2.sif"
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
