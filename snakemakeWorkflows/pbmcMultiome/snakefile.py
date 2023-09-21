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


rule singlecells_wnn_analysis_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
        packages=  "config/installedAtacPackages"
  output: "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
          "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
          "output/pbmcMultiome/singlecells_analysis/seurat.SCT.h5ad",
          "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad",
  conda: SuperCellMultiomics
  shell: "Rscript R/wnnAnalysis10xMultiomeCL.R -i pbmcMultiome -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis \
          -p 1:40 -q 2:40 -a SCTransform"

rule metacell_identification_pbmc_multiome:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellMulti/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:40 -q 2:40 -r SCT -e TRUE -f TRUE -k 30 -g {wildcards.gamma}"
          
rule metacell_identification_pbmc_RNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellRNA/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  benchmark : "benchmark/pbmcMultiome/SuperCellRNA/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellRNA/g{wildcards.gamma}/ \
          -p 1:40 -r SCT -e TRUE -f TRUE -k 30 -g {wildcards.gamma}"
          

          
rule metacell_identification_pbmc_ATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellATAC/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  benchmark : "benchmark/pbmcMultiome/SuperCellATAC/g{gamma}/seurat.multiome.mc.txt"
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellATAC/g{wildcards.gamma}/ \
          -q 2:40 -e TRUE -f TRUE -k 30 -g {wildcards.gamma}"

rule metacell_identification_pbmc_seacellsRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.SCT.h5ad"
  output: "output/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  benchmark:"benchmark/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsRNA/g{wildcards.gamma}/ \
          -d 1:40 -r pca -g {wildcards.gamma}"
          
rule data_aggregation_pbmc_seacellsRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsRNA/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsRNA/g{wildcards.gamma}/ \
          -f TRUE -g {wildcards.gamma} -x SEACell-"
          
rule metacell_identification_pbmc_MetaCellRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  benchmark: "benchmark/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.txt"
  conda: SuperCellMultiomics_pyenv
  shell: "python3 python/MetaCell2CL.py -i {input.singlecells} \
          -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"
          
rule data_aggregation_pbmc_MetaCellRNA:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
          -f TRUE -g {wildcards.gamma}"
          
          
rule metacell_identification_pbmc_seacellsATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
  output: "output/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  benchmark:"benchmark/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsATAC/g{wildcards.gamma}/ \
          -d 2:40 -r lsi -g {wildcards.gamma}"
          
rule data_aggregation_pbmc_seacellsATAC:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsATAC/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsATAC/g{wildcards.gamma}/ \
          -f TRUE -g {wildcards.gamma} -x SEACell-"
          
rule random_metacell_pbmc_multiome:
  input: 
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/randomMetacells/g{gamma}/seurat.multiome.mc.rds"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiomeCL.R -i {input.singlecells} -d TRUE \
          -o output/pbmcMultiome/randomMetacells/g{wildcards.gamma}/ \
          -f TRUE -g {wildcards.gamma}"

          
rule gene_selection_for_correlation:
  input: "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/selectedGenes.txt"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'seurat.sc <- readRDS(\"{input}\");\
          selectedGenes <- rownames(seurat.sc[[\"RNA\"]]@counts)[rowSums(seurat.sc[[\"RNA\"]]@counts)> ncol(seurat.sc)*0.005];\
          write.table(selectedGenes,\"{output}\")'"

rule benchmark_and_correlations:
  input: metacells = "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mc.rds",
         geneList = "output/pbmcMultiome/selectedGenes.txt",
         singleCells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTable.csv",
          "output/pbmcMultiome/{inputMetacells}/g{gamma}/chromVarCorrTable.csv",
          "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mc.activities.rds"
  params: outdir = "output/pbmcMultiome/{inputMetacells}/g{gamma}/",
          python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript R/atacRnaCorrAnalysisCL.R -i {input.metacells} \
  -g {input.geneList} \
  -l {input.singleCells} \
  -o {params.outdir} \
  -p TRUE -e {params.python} \
  -s TRUE -c TRUE -w 10 \
  -n 1:40 -q 2:40"

rule atacRnaCorrSingleCellAnalysis:
  input: metacells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
         geneList = "output/pbmcMultiome/selectedGenes.txt"
  output: "output/pbmcMultiome/singlecells_analysis/corrTable.csv",
          "output/pbmcMultiome/singlecells_analysis/chromVarCorrTable.csv",
          "output/pbmcMultiome/singlecells_analysis/seurat.multiome.mc.activities.rds"
  params: outdir = "output/pbmcMultiome/singlecells_analysis"
  conda: SuperCellMultiomics
  shell: "Rscript R/atacRnaCorrAnalysisCL.R -i {input.metacells}  -g {input.geneList} -o {params.outdir} -s FALSE -c TRUE"
