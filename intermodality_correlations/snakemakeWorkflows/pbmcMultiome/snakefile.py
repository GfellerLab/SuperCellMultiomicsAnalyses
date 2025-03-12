#################################################################################################
#################### PBMC 10X multiome dataset ##################################################
#################################################################################################

fragmentFile = config["pbmcMultiome"]["fragmentFile"]


rule download_fragment_file_pbmc_multiome:
  output: "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
           "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
  params: file_url = fragmentFile
  shell: "cd input/pbmcMultiome/;\
          wget {params.file_url};\
          wget {params.file_url}.tbi"
          
          
rule install_seurat_pbmcmultiome:
  input : frags =  "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
  output: "input/pbmcMultiome/pbmcMultiome_broad_peaks.rds"
  singularity: "config/matk_multiomics.sif"
  shell: "Rscript -e 'library(Seurat);\
                      library(Seurat);library(Signac);\
                      library(dplyr);library(EnsDb.Hsapiens.v86);\
                      library(BSgenome.Hsapiens.UCSC.hg38);\
                      .libPaths(c(\"./input/\", .libPaths()));\
                      SeuratData::InstallData(\"pbmcMultiome\",lib = \"./input/\");\
                      library(SeuratData); \
                      data(\"pbmc.atac\"); \
                      grange.counts <- StringToGRanges(rownames(pbmc.atac)) ;\
                      grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts); \
                      pbmc.atac <- pbmc.atac[as.vector(grange.use), ]; \
                      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86); \
                      seqlevelsStyle(annotations) <- \"UCSC\"; \
                      genome(annotations) <- \"hg38\"; \
                      message(\"creating chromatine assay\"); \
                      pbmc.atac <- CreateChromatinAssay(counts = pbmc.atac@assays$ATAC@counts, \
                                                        genome = \"hg38\", \
                                                        fragments = \"{input.frags}\", \
                                                        annotation = annotations \
                                                        ); \
                      data(\"pbmc.rna\");  \
                      pbmc.rna[[\"ATAC\"]] <- pbmc.atac;  \
                      remove(pbmc.atac);  \
                      pbmc <- pbmc.rna;  \
                      remove(pbmc.rna);  \
                      pbmc.rna <- data(\"pbmc.rna\"); \
                      saveRDS(pbmc,\"{output}\")'"
                      
rule singlecells_macs2_peak_calling_pbmc_multiome:
  input: "input/pbmcMultiome/pbmcMultiome_broad_peaks.rds"
  output: "input/pbmcMultiome/pbmcMultiome.rds"
  singularity: "config/matk_multiomics.sif"
  shell: """
         Rscript snakemakeWorkflows/pbmcMultiome/R/macs2PeakCallingCL.R \
          -i {input} \
          -o input/pbmcMultiome/
         """
                      
          
rule singlecells_rna_pca_analysis_pbmc_multiome:
  input:  dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  singularity: "config/matk_multiomics.sif"
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seurat.RNA.txt"
  shell: "Rscript R/RNA_pca_analysis10xMultiomeCL.R -i {input.dataset}  \
          -o output/pbmcMultiome/singlecells_analysis \
          -a SCTransform"
          
rule singlecells_atac_lsi_analysis_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
         index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
         dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
  singularity: "config/matk_multiomics.sif"
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seurat.ATAC.txt"
  shell: "Rscript R/ATAC_lsi_analysis10xMultiomeCL.R -i {input.dataset}\
          -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis"

rule preprocessing_supercell_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
        dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds",
  singularity: "config/matk_multiomics.sif"
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seurat_multimodal.txt"
  shell: "Rscript R/preprocessing_seurat_for_SuperCell_10xMultiomeCL.R -i {input.dataset} -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis \
          -a SCTransform"
          
rule wnn_single_cell_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
        dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
  singularity: "config/matk_multiomics.sif"
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seuratWNN.txt"
  shell: "Rscript R/wnnAnalysis10xMultiomeCL.R -i {input.dataset} -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis \
          -d -y {params.python} \
          -p 1:40 -q 2:40 -a SCTransform "




rule metacell_identification_pbmc_multiome:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds"
  output: "output/pbmcMultiome/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:40 -q 2:40 -r SCT -e TRUE -k 30 -g {wildcards.gamma} -b"

rule data_aggregation_pbmc_SuperCellMulti:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  output: "output/pbmcMultiome/SuperCellMulti/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/SuperCellMulti/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"
          
rule metacell_identification_pbmc_RNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/SuperCellRNA/g{gamma}/SuperCellMemberships.csv"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  benchmark : "benchmark/pbmcMultiome/SuperCellRNA/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellRNA/g{wildcards.gamma}/ \
          -p 1:40 -r SCT -e TRUE -k 30 -g {wildcards.gamma} -b"
          
rule data_aggregation_pbmc_SuperCellRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/SuperCellRNA/g{gamma}/SuperCellMemberships.csv"
  output: "output/pbmcMultiome/SuperCellRNA/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/SuperCellRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"
          
rule metacell_identification_pbmc_ATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad"
  output: "output/pbmcMultiome/SuperCellATAC/g{gamma}/SuperCellMemberships.csv"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  benchmark : "benchmark/pbmcMultiome/SuperCellATAC/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellATAC/g{wildcards.gamma}/ \
          -q 2:40 -e TRUE -k 30 -g {wildcards.gamma} -b"

rule data_aggregation_pbmc_SuperCellATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/SuperCellATAC/g{gamma}/SuperCellMemberships.csv"
  output: "output/pbmcMultiome/SuperCellATAC/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/SuperCellATAC/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

rule metacell_identification_pbmc_seacellsRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.csv"
  singularity: "config/matk_multiomics.sif"
  benchmark:"benchmark/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsRNA/g{wildcards.gamma}/ \
          -d 1:40 -r pca -g {wildcards.gamma}"
          
rule data_aggregation_pbmc_seacellsRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsRNA/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x SEACell-"
          
# rule metacell_identification_pbmc_MetaCellRNA:
  # input:
  #       singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  # output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  # benchmark: "benchmark/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.txt"
  # conda: SuperCellMultiomics_pyenv
#   shell: "python3 python/MetaCell2CL.py -i {input.singlecells} \
#           -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
#           -g {wildcards.gamma}"
          
rule metacell_identification_pbmc_MetaCellRNA:
  input: 
                singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  benchmark: "benchmark/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.txt"
  singularity: "config/matk_multiomics.sif"
  shell: "python3 python/MATK_MetaCell2CL.py -i {input.singlecells} \
          -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"
          
rule data_aggregation_pbmc_MetaCellRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"
          
          
rule metacell_identification_pbmc_seacellsATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
  output: "output/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.csv"
  singularity: "config/matk_multiomics.sif"
  benchmark:"benchmark/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsATAC/g{wildcards.gamma}/ \
          -d 2:40 -r lsi -g {wildcards.gamma}"
          
rule data_aggregation_pbmc_seacellsATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsATAC/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsATAC/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x SEACell-"
          
rule random_metacell_pbmc_multiome:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/randomMetacells/g{gamma}/seurat.multiome.mc.rds"
  singularity: "config/matk_multiomics.sif"
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -d TRUE \
          -o output/pbmcMultiome/randomMetacells/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

          
rule gene_selection_for_correlation:
  input: "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/selectedGenes.txt"
  singularity: "config/matk_multiomics.sif"
  shell: "Rscript -e 'seurat.sc <- readRDS(\"{input}\");\
          selectedGenes <- rownames(seurat.sc[[\"RNA\"]]@counts)[rowSums(seurat.sc[[\"RNA\"]]@counts)> ncol(seurat.sc)*0.005];\
          write.table(selectedGenes,\"{output}\")'"
          
  
        
####################################################################################################################################
########################################## Benchmark and Correlations  #############################################################
####################################################################################################################################



rule benchmark_and_correlations:
  input: metacells = "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mc.rds",
         geneList = "output/pbmcMultiome/selectedGenes.txt",
         singleCells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
         motifs =  "input/HOCCOMOCO_motifs/H12CORE_pfm_annotated.rds"
  output: "output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTable.csv",
          "output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTablePearson.csv",
          "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mc.activities.rds"
  params: outdir = "output/pbmcMultiome/{inputMetacells}/g{gamma}/",
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  singularity: "config/matk_multiomics.sif"
  shell: "Rscript R/atacRnaCorrAnalysis_v5_CL.R \
  -i {input.metacells} \
  -g {input.geneList} \
  -l {input.singleCells} \
  -o {params.outdir} \
  -c {input.motifs} \
  -p -e {params.python} \
  -s TRUE -w 10 \
  -n 1:40 -q 2:40"

rule atacRnaCorrSingleCellAnalysis:
  input: metacells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
         geneList = "output/pbmcMultiome/selectedGenes.txt",
         motifs =  "input/HOCCOMOCO_motifs/H12CORE_pfm_annotated.rds"
  output: "output/pbmcMultiome/singlecells_analysis/corrTable.csv",
          "output/pbmcMultiome/singlecells_analysis/corrTablePearson.csv",
          "output/pbmcMultiome/singlecells_analysis/seurat.multiome.mc.activities.rds"
  params: outdir = "output/pbmcMultiome/singlecells_analysis",
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  singularity: "config/matk_multiomics.sif"
  shell: "Rscript R/atacRnaCorrAnalysis_v5_CL.R -i {input.metacells}  -g {input.geneList} -o {params.outdir} -c {input.motifs} -s FALSE"

rule report_all_bench_corr:
  input: "reports/pbmcMultiome/pbmcMultiome.Rmd",
        "output/pbmcMultiome/singlecells_analysis/seurat.multiome.mc.activities.rds",
        expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTablePearson.csv",gamma = GAMMA,inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"])
  output: "reports/pbmcMultiome/pbmcMultiome.html"
  singularity: "config/matk_multiomics.sif"
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"

  
  
  
