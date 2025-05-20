#############################################
##########     Kernel tests      ############
#############################################

  
rule download_mixology:
  output: 'input/sc_mixology/data/sincell_with_class_5cl.RData'
  shell: "cd input; git clone https://github.com/LuyiTian/sc_mixology.git"
  
rule load_mixology:
  input: "input/sc_mixology/data/sincell_with_class_5cl.RData"
  output: "input/sce_sc_10x_5cl_qc.rds"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'load(\"{input}\");\
          sce_sc_10x_5cl_qc <- Seurat::as.Seurat(sce_sc_10x_5cl_qc);\
          sce_sc_10x_5cl_qc <- Seurat::NormalizeData(sce_sc_10x_5cl_qc);\
          sce_sc_10x_5cl_qc <- SeuratObject::RenameAssays(object = sce_sc_10x_5cl_qc, originalexp = \"RNA\");\
          saveRDS(sce_sc_10x_5cl_qc,\"{output}\")'"
          
          
rule rarePopDetectionGamma_mixology:
  input: 'input/sce_sc_10x_5cl_qc.rds'
  output: 'output/testKernel/rarePop_mixology/{cellType}/detectionGammaRes.csv'
  params: outputDir = "output/testKernel/rarePop_mixology/{cellType}/"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionGammaCL.R -i {input} -o {params.outputDir} -c cell_line -t {wildcards.cellType} -r 0.005 -p 1:30 -n 30 -v 2000"
  
rule rarePopDetectionSize_mixology:
  wildcard_constraints: cellType= '\w+',
  input: 'input/sce_sc_10x_5cl_qc.rds'
  output: 'output/testKernel/rarePop_mixology/{cellType}/kernel_{kernel}/detectionSizeRes.csv'
  params: outputDir = "output/testKernel/rarePop_mixology/{cellType}/kernel_{kernel}"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionSizeCL.R -i {input} -o {params.outputDir} -c cell_line -t {wildcards.cellType} -r 0.005 -p 1:30 -n 30 -v 2000 -d {wildcards.kernel} -g 20"


rule rarePopDetectionSize_mixology_k5:
  input: 'input/sce_sc_10x_5cl_qc.rds'
  output: 'output/testKernel/rarePop_mixology/{cellType}/kernel_FALSE_k5/detectionSizeRes.csv'
  params: outputDir = "output/testKernel/rarePop_mixology/{cellType}/kernel_FALSE_k5"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionSizeCL.R -i {input} -o {params.outputDir} -c cell_line -t {wildcards.cellType} -r 0.005 -p 1:30 -n 30 -v 2000 -d FALSE -g 20 -k 5 -w 8"


rule rarePopDetectionGamma_bm_cite_seq:
  input:  "config/installed_bmcite"
  output: "output/testKernel/CITEseq/bmcite/{BMcellType}/detectionGammaRes.csv"
  params: outputDir = "output/testKernel/CITEseq/bmcite/{BMcellType}/"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionGammaCiteSeqCL.R -i bmcite -o {params.outputDir} -c celltype.l2 -t {wildcards.BMcellType} -p 1:30 -q 1:18 -v 2000 -r 0.001 -n 30"
  
  
rule rarePopDetectionSize_bm_cite_seq:
  input:  "config/installed_bmcite"
  output: "output/testKernel/CITEseq/bmcite/{BMcellType}/kernel_{kernel}/detectionSizeRes.csv"
  params: outputDir = "output/testKernel/CITEseq/bmcite/{BMcellType}/kernel_{kernel}"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionSizeCiteSeqCL.R -i bmcite -o {params.outputDir} -c celltype.l2 -t {wildcards.BMcellType} -p 1:30 -q 1:18 -n 30 -v 2000 -d {wildcards.kernel} -g 20 -w 20"



rule rarePopDetectionGamma_pbmc_multiome:
  input:  "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/testKernel/10Xmultiome/PBMC/{PBMCcellType}/detectionGammaRes.csv"
  params: outputDir = "output/testKernel/10Xmultiome/PBMC/{PBMCcellType}/"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionGammaMultiomeCL.R -i {input}\
          -a SCT\
          -o {params.outputDir} -c seurat_annotations -t {wildcards.PBMCcellType} -p 1:40 -q 2:40 -r 0.001 -n 30"
  
  
rule rarePopDetectionSize_pbmc_multiome:
  input:  "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/testKernel/10Xmultiome/PBMC/{PBMCcellType}/kernel_{kernel}/detectionSizeRes.csv"
  params: outputDir = "output/testKernel/10Xmultiome/PBMC/{PBMCcellType}/kernel_{kernel}"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/testKernel/R/rarePopDetectionSizeMultiomeCL.R -i {input}\
         -a SCT\
         -o {params.outputDir} -c seurat_annotations -t {wildcards.PBMCcellType} -p 1:40 -q 2:40 -n 30 -v 2000 -d {wildcards.kernel} -g 20 -w 6"

rule detect_pop_multi:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds"
  output: "output/testKernel/pbmcMultiome/detectRarePop/g{gamma}/seed_{seed}.csv"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript snakemakeWorkflows/testKernel/R/detect_rare_pop_multi_for_loop.R -i {input.singlecells} \
          -o output/testKernel/pbmcMultiome/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -s {wildcards.seed}"

rule all_seed_detect_rare_pop_multiome:
  input: expand("output/testKernel/pbmcMultiome/detectRarePop/g{gamma}/seed_{seed}.csv",gamma = [20,75], seed = range(1,30,1))
  output: "output/testKernel/pbmcMultiome/detectRarePop/all_seeds_processed_at_all_gammas"
  shell: "touch {output}"

          
          
rule reportRarePopExp:
  input: "reports/testKernel/kernel_analyses.Rmd",
         expand("config/installed_{seuratDataset}",seuratDataset = seuratDatasets),
         expand("output/testKernel/CITEseq/bmcite/{BMcellType}/detectionGammaRes.csv",BMcellType = cellTypesBmCiteSeq),
         expand("output/testKernel/CITEseq/bmcite/{BMcellType}/kernel_{kernel}/detectionSizeRes.csv",BMcellType = cellTypesBmCiteSeq,kernel = ["FALSE","TRUE"]),
         expand("output/testKernel/10Xmultiome/PBMC/{PBMCcellType}/detectionGammaRes.csv",PBMCcellType = cellTypesPBMCcellType10Xmultiome),
         expand("output/testKernel/10Xmultiome/PBMC/{PBMCcellType}/kernel_{kernel}/detectionSizeRes.csv",PBMCcellType = cellTypesPBMCcellType10Xmultiome,kernel = ["FALSE","TRUE"]),
         expand("output/testKernel/rarePop_mixology/{cellType}/detectionGammaRes.csv",cellType = cellTypesMixology),
         expand("output/testKernel/rarePop_mixology/{cellType}/kernel_{kernel}/detectionSizeRes.csv",cellType = cellTypesMixology,kernel = ["FALSE","TRUE"]),
         expand("output/testKernel/rarePop_mixology/{cellType}/kernel_FALSE_k5/detectionSizeRes.csv",cellType = cellTypesMixology)
  output: "reports/testKernel/kernel_analyses.html"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
