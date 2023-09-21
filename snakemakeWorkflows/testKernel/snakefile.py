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


rule reportRarePopExp:
  input: "reports/testKernel/kernel_analyses.Rmd",
         expand("config/installed_{seuratDataset}",seuratDataset = seuratDatasets),
         expand("output/testKernel/CITEseq/bmcite/{BMcellType}/detectionGammaRes.csv",BMcellType = cellTypesBmCiteSeq),
         expand("output/testKernel/CITEseq/bmcite/{BMcellType}/kernel_{kernel}/detectionSizeRes.csv",BMcellType = cellTypesBmCiteSeq,kernel = ["FALSE","TRUE"]),
         expand("output/testKernel/rarePop_mixology/{cellType}/detectionGammaRes.csv",cellType = cellTypesMixology),
         expand("output/testKernel/rarePop_mixology/{cellType}/kernel_{kernel}/detectionSizeRes.csv",cellType = cellTypesMixology,kernel = ["FALSE","TRUE"]),
         expand("output/testKernel/rarePop_mixology/{cellType}/kernel_FALSE_k5/detectionSizeRes.csv",cellType = cellTypesMixology)
  output: "reports/testKernel/kernel_analyses.html"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
