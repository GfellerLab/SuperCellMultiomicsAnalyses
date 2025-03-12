#################################################################################################
#################### hspc 10X multiome dataset ##################################################
#################################################################################################

fragmentFileRep1 = config["hspcMultiomePersad"]["fragmentFileRep1"]
fragmentFileRep2 = config["hspcMultiomePersad"]["fragmentfileRep2"]
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
  input:  packages=  "config/installedAtacPackages",
          adata = "input/hspcMultiomePersad/cd34_multiome_rna.h5ad"
  output: "output/hspcMultiomePersad/singlecells_analysis/seurat.RNA.h5ad"
  conda: SuperCellMultiomics
  benchmark: "benchmark/hspcMultiomePersad/singlecells_analysis/seurat.RNA.txt"
  shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/RNA_pca_hspcMultiomeCL.R -i {input.adata} \
          -o output/hspcMultiomePersad/singlecells_analysis"
          
rule singlecells_atac_lsi_analysis_hspc_multiome:
  input: fragment_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz",
         fragment_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz",
         index_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz.tbi",
         index_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz.tbi",
         adata = "input/hspcMultiomePersad/cd34_multiome_atac.h5ad",
         packages=  "config/installedAtacPackages"
  output: "output/hspcMultiomePersad/singlecells_analysis/seurat.ATAC.h5ad",
  conda: SuperCellMultiomics
  benchmark: "benchmark/hspcMultiomePersad/singlecells_analysis/seurat.ATAC.txt"
  shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/ATAC_lsi_hspcMultiomeCL.R -i {input.adata}\
          -f {input.fragment_rep1}+{input.fragment_rep2} \
          -o output/hspcMultiomePersad/singlecells_analysis"

rule preprocessing_supercell_hspc_multiome:
  input: fragment_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz",
         fragment_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz",
         index_rep1 = "input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz.tbi",
         index_rep2 = "input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz.tbi",
         atac = "input/hspcMultiomePersad/cd34_multiome_atac.h5ad",
         rna = "input/hspcMultiomePersad/cd34_multiome_rna.h5ad",
         packages=  "config/installedAtacPackages"
  output: "output/hspcMultiomePersad/singlecells_analysis/seurat_multimodal.rds",
  conda: SuperCellMultiomics
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
         rna = "input/hspcMultiomePersad/cd34_multiome_rna.h5ad",
         packages=  "config/installedAtacPackages"
  output: "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds","output/hspcMultiomePersad/singlecells_analysis/selectedGenes.txt"
  conda: SuperCellMultiomics
  benchmark: "benchmark/hspcMultiomePersad/singlecells_analysis/seuratWNN.txt"
  shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/wnnAnalysis10xMultiomeCL.R -i {input.rna}\
          -j {input.atac} -f {input.fragment_rep1}+{input.fragment_rep2} \
          -p 1:50 -q 2:50 \
          -o output/hspcMultiomePersad/singlecells_analysis"

# rule singlecells_macs2_peak_calling_per_cluster_hspc_multiome:
#   input: "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds"
#   output: "output/hspcMultiomePersad/singlecells_analysis/seuratWNN_with_macs2_peaks.rds"
#   conda: SuperCellMultiomics
#   shell: "Rscript snakemakeWorkflows/hspcMultiomePersad/R/macs2PeakCallingPerSingleCellClustersCL.R\
#           -i {input}\
#           -o output/hspcMultiomePersad/singlecells_analysis/"


rule metacell_identification_hspc_multiome:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat_multimodal.rds"
  output: "output/hspcMultiomePersad/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  conda: SuperCellMultiomics
  params: workdir = wdir
  benchmark :  "benchmark/hspcMultiomePersad/SuperCellMulti/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} \
          -o output/hspcMultiomePersad/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:50 -q 2:50 -e TRUE -k 30 -g {wildcards.gamma} -b"

rule data_aggregation_hspc_SuperCellMulti:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  output: "output/hspcMultiomePersad/SuperCellMulti/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/SuperCellMulti/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"
          
rule metacell_identification_hspc_RNA:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/hspcMultiomePersad/SuperCellRNA/g{gamma}/SuperCellMemberships.csv"
  conda: SuperCellMultiomics
  params: workdir = wdir
  benchmark : "benchmark/hspcMultiomePersad/SuperCellRNA/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} \
          -o output/hspcMultiomePersad/SuperCellRNA/g{wildcards.gamma}/ \
          -p 1:50 -e TRUE -k 30 -g {wildcards.gamma} -b"
          
rule data_aggregation_hspc_SuperCellRNA:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/SuperCellRNA/g{gamma}/SuperCellMemberships.csv"
  output: "output/hspcMultiomePersad/SuperCellRNA/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/SuperCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"
          
rule metacell_identification_hspc_ATAC:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.ATAC.h5ad"
  output: "output/hspcMultiomePersad/SuperCellATAC/g{gamma}/SuperCellMemberships.csv"
  conda: SuperCellMultiomics
  params: workdir = wdir
  benchmark : "benchmark/hspcMultiomePersad/SuperCellATAC/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} \
          -o output/hspcMultiomePersad/SuperCellATAC/g{wildcards.gamma}/ \
          -q 2:50 -e TRUE -k 30 -g {wildcards.gamma} -b"

rule data_aggregation_hspc_SuperCellATAC:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/SuperCellATAC/g{gamma}/SuperCellMemberships.csv"
  output: "output/hspcMultiomePersad/SuperCellATAC/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/SuperCellATAC/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"

rule metacell_identification_hspc_seacellsRNA:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/hspcMultiomePersad/seacellsRNA/g{gamma}/seacellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  benchmark:"benchmark/hspcMultiomePersad/seacellsRNA/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/hspcMultiomePersad/seacellsRNA/g{wildcards.gamma}/ \
          -d 1:50 -r pca -g {wildcards.gamma}"
          
rule data_aggregation_hspc_seacellsRNA:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/seacellsRNA/g{gamma}/seacellMemberships.csv"
  output: "output/hspcMultiomePersad/seacellsRNA/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/seacellsRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -x SEACell-"
          
rule metacell_identification_hspc_MetaCellRNA:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/hspcMultiomePersad/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  benchmark: "benchmark/hspcMultiomePersad/MetaCellRNA/g{gamma}/MetaCellMemberships.txt"
  conda: metacell2_0_9_env
  shell: "python3 python/MATK_MetaCell2CL.py -i {input.singlecells} \
          -o output/hspcMultiomePersad/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"
          
rule data_aggregation_hspc_MetaCellRNA:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  output: "output/hspcMultiomePersad/MetaCellRNA/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"
          
          
rule metacell_identification_hspc_seacellsATAC:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seurat.ATAC.h5ad",
  output: "output/hspcMultiomePersad/seacellsATAC/g{gamma}/seacellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  benchmark:"benchmark/hspcMultiomePersad/seacellsATAC/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/hspcMultiomePersad/seacellsATAC/g{wildcards.gamma}/ \
          -d 2:50 -r lsi -g {wildcards.gamma}"
          
rule data_aggregation_hspc_seacellsATAC:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
        memberships = "output/hspcMultiomePersad/seacellsATAC/g{gamma}/seacellMemberships.csv"
  output: "output/hspcMultiomePersad/seacellsATAC/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/hspcMultiomePersad/seacellsATAC/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -x SEACell-"
          
rule random_metacell_hspc_multiome:
  input: 
        singlecells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds"
  output: "output/hspcMultiomePersad/randomMetacells/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -d TRUE \
          -o output/hspcMultiomePersad/randomMetacells/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"

          
# rule gene_selection_for_correlation_hspc:
#   input: "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds"
#   output: 
#   conda: SuperCellMultiomics
#   shell: "Rscript -e 'seurat.sc <- readRDS(\"{input}\");\
#           selectedGenes <- rownames(seurat.sc[[\"RNA\"]]@counts)[Matrix::rowSums(seurat.sc[[\"RNA\"]]@counts)> ncol(seurat.sc)*0.005];\
#           write.table(selectedGenes,\"{output}\")'"
          

        
        
####################################################################################################################################
########################################## Benchmark and Correlations  #############################################################
####################################################################################################################################



rule benchmark_and_correlations_hspc:
  input: metacells = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.mc.rds",
         geneList = "output/hspcMultiomePersad/singlecells_analysis/selectedGenes.txt",
         singleCells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
         motifs =  "input/motifs/human/jaspar2024_pfm.rds"
  output: "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/corrTable.csv",
          "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/corrTablePearson.csv",
          "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/chromVarCorrTable.csv",
          "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/seurat.multiome.mc.activities.rds"
  params: outdir = "output/hspcMultiomePersad/{inputMetacells}/g{gamma}/",
          python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript R/atacRnaCorrAnalysisCL.R \
  -i {input.metacells} \
  -g {input.geneList} \
  -l {input.singleCells} \
  -o {params.outdir} \
  -c {input.motifs} \
  -p -e {params.python} \
  -s TRUE -w 10 \
  -n 1:50 -q 2:50"

rule atacRnaCorrSingleCellAnalysis_hspc:
  input: metacells = "output/hspcMultiomePersad/singlecells_analysis/seuratWNN.rds",
         geneList = "output/hspcMultiomePersad/singlecells_analysis/selectedGenes.txt",
         motifs =  "input/motifs/human/jaspar2024_pfm.rds"
  output: "output/hspcMultiomePersad/singlecells_analysis/corrTable.csv",
          "output/hspcMultiomePersad/singlecells_analysis/corrTablePearson.csv",
          "output/hspcMultiomePersad/singlecells_analysis/chromVarCorrTable.csv",
          "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.mc.activities.rds"
  params: outdir = "output/hspcMultiomePersad/singlecells_analysis"
  conda: SuperCellMultiomics
  shell: "Rscript R/atacRnaCorrAnalysisCL.R -i {input.metacells}  -g {input.geneList} -o {params.outdir} -c {input.motifs} -s FALSE"

rule report_all_bench_corr_hspc:
  input: "reports/hspcMultiomePersad/hspcMultiomePersad.Rmd",
         "output/hspcMultiomePersad/singlecells_analysis/seurat.multiome.mc.activities.rds",
         expand("output/hspcMultiomePersad/{inputMetacells}/g{gamma}/corrTablePearson.csv",gamma = GAMMA,inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"])
  output: "reports/hspcMultiomePersad/hspcMultiomePersad.html"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"

  
  
  
