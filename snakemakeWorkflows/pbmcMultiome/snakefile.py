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
  singularity: config["sif_file"]
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
  singularity: config["sif_file"]
  shell: """
         Rscript snakemakeWorkflows/pbmcMultiome/R/macs2PeakCallingCL.R \
          -i {input} \
          -o input/pbmcMultiome/
         """

rule singlecells_rna_pca_analysis_pbmc_multiome:
  input:  dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  singularity: config["sif_file"]
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seurat.RNA.txt"
  shell: "Rscript R/RNA_pca_analysis10xMultiomeCL.R -i {input.dataset}  \
          -o output/pbmcMultiome/singlecells_analysis \
          -a SCTransform"

rule singlecells_atac_lsi_analysis_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
         index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
         dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
  singularity: config["sif_file"]
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seurat.ATAC.txt"
  shell: "Rscript R/ATAC_lsi_analysis10xMultiomeCL.R -i {input.dataset}\
          -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis"

rule preprocessing_supercell_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
        dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds",
  singularity: config["sif_file"]
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seurat_multimodal.txt"
  shell: "Rscript R/preprocessing_seurat_for_SuperCell_10xMultiomeCL.R -i {input.dataset} -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis \
          -a SCTransform"

rule wnn_single_cell_pbmc_multiome:
  input: fragment = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        index = "input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
        dataset=  "input/pbmcMultiome/pbmcMultiome.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
  singularity: config["sif_file"]
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/seuratWNN.txt"
  shell: "Rscript R/wnnAnalysis10xMultiomeCL.R -i {input.dataset} -f {input.fragment} \
          -o output/pbmcMultiome/singlecells_analysis \
          -d -y {params.python} \
          -p 1:40 -q 2:40 -a SCTransform "

rule singlecells_mofa_analysis_pbmc_multiome:
  input: dataset=  "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
         script= "R/mofaAnalysis10xMultiomeCL.R"
  output: "output/pbmcMultiome/singlecells_analysis/adata_mofa.h5ad"
  singularity: config["sif_file_mofa"]
  params: python = "/opt/conda/envs/mofa_env/bin/python"
  benchmark: "benchmark/pbmcMultiome/singlecells_analysis/adata_mofa.txt"
  shell: "Rscript {input.script} \
          -i {input.dataset} \
          -o output/pbmcMultiome/singlecells_analysis \
          -n 40 \
          -y {params.python}"

rule metacell_identification_pbmc_multiome:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds"
  output: "output/pbmcMultiome/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
  singularity: config["sif_file"]
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
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/SuperCellMulti/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

rule metacell_identification_pbmc_RNA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/SuperCellRNA/g{gamma}/SuperCellMemberships.csv"
  singularity: config["sif_file"]
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
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/SuperCellRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

rule metacell_identification_pbmc_ATAC:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad"
  output: "output/pbmcMultiome/SuperCellATAC/g{gamma}/SuperCellMemberships.csv"
  singularity: config["sif_file"]
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
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/SuperCellATAC/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"
          
rule metacell_identification_pbmc_seacellsMOFA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/adata_mofa.h5ad"
  output: "output/pbmcMultiome/seacellsMOFA/g{gamma}/seacellMemberships.csv"
  singularity: config["sif_file"]
  benchmark:"benchmark/pbmcMultiome/seacellsMOFA/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsMOFA/g{wildcards.gamma}/ \
          -r mofa -g {wildcards.gamma}"

rule data_aggregation_pbmc_seacellsMOFA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsMOFA/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsMOFA/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsMOFA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x SEACell-"
          
rule metacell_identification_pbmc_metaqRNA:
  input:
      rna = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/metaqRNA/g{gamma}/metaqMemberships.csv"
  singularity: config["sif_file_metaq"]
  benchmark: "benchmark/pbmcMultiome/metaqRNA/g{gamma}/benchIdentification.txt"
  shell: "/opt/conda/envs/MetaQ/bin/python python/MetaQ_CL.py --data_path {input.rna}\
         --data_type 'RNA' \
         --outdir output/pbmcMultiome/metaqRNA/g{wildcards.gamma}/ \
         --gamma {wildcards.gamma}"

rule data_aggregation_pbmc_metaqRNA:
  input:
       singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
       memberships = "output/pbmcMultiome/metaqRNA/g{gamma}/metaqMemberships.csv"
  output: "output/pbmcMultiome/metaqRNA/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"] 
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/metaqRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x MetaQ-"
          
rule metacell_identification_pbmc_metaqATAC:
  input:
      rna = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad"
  output: "output/pbmcMultiome/metaqATAC/g{gamma}/metaqMemberships.csv"
  singularity: config["sif_file_metaq"]
  benchmark: "benchmark/pbmcMultiome/metaqATAC/g{gamma}/benchIdentification.txt"
  shell: "/opt/conda/envs/MetaQ/bin/python python/MetaQ_CL.py --data_path {input.rna}\
         --data_type 'RNA' \
         --outdir output/pbmcMultiome/metaqATAC/g{wildcards.gamma}/ \
         --gamma {wildcards.gamma}"

rule data_aggregation_pbmc_metaqATAC:
  input:
       singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
       memberships = "output/pbmcMultiome/metaqATAC/g{gamma}/metaqMemberships.csv"
  output: "output/pbmcMultiome/metaqATAC/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"] 
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/metaqATAC/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x MetaQ-"
          

rule metacell_identification_pbmc_metaqMulti:
  input:
      atac = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
      rna = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/metaqMulti/g{gamma}/metaqMemberships.csv"
  singularity: config["sif_file_metaq"]
  benchmark: "benchmark/pbmcMultiome/metaqMulti/g{gamma}/benchIdentification.txt"
  shell: "/opt/conda/envs/MetaQ/bin/python python/MetaQ_CL.py --data_path {input.rna} {input.atac} \
         --data_type 'RNA' 'ATAC' \
         --outdir output/pbmcMultiome/metaqMulti/g{wildcards.gamma}/ \
         --gamma {wildcards.gamma}"

rule data_aggregation_pbmc_metaqMulti:
  input:
       singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
       memberships = "output/pbmcMultiome/metaqMulti/g{gamma}/metaqMemberships.csv"
  output: "output/pbmcMultiome/metaqMulti/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"] 
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/metaqMulti/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x MetaQ-"

rule metacell_identification_pbmc_seacellsRNA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.csv"
  singularity: config["sif_file"]
  benchmark:"benchmark/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsRNA/g{wildcards.gamma}/ \
          -d 1:40 -r pca -g {wildcards.gamma}"

rule data_aggregation_pbmc_seacellsRNA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsRNA/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsRNA/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x SEACell-"


rule metacell_identification_pbmc_MetaCellRNA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  benchmark: "benchmark/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.txt"
  singularity: config["sif_file"]
  shell: "python3 python/MATK_MetaCell2CL.py -i {input.singlecells} \
          -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma}"

rule data_aggregation_pbmc_MetaCellRNA:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  output: "output/pbmcMultiome/MetaCellRNA/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/MetaCellRNA/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"


rule metacell_identification_pbmc_seacellsATAC:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.ATAC.h5ad",
  output: "output/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.csv"
  singularity: config["sif_file"]
  benchmark:"benchmark/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcMultiome/seacellsATAC/g{wildcards.gamma}/ \
          -d 2:40 -r lsi -g {wildcards.gamma}"

rule data_aggregation_pbmc_seacellsATAC:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
        memberships = "output/pbmcMultiome/seacellsATAC/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcMultiome/seacellsATAC/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcMultiome/seacellsATAC/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma} -x SEACell-"

rule random_metacell_pbmc_multiome:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/randomMetacells/g{gamma}/seurat.multiome.mc.rds"
  singularity: config["sif_file"]
  params: workdir = wdir
  shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
          Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -d TRUE \
          -o output/pbmcMultiome/randomMetacells/g{wildcards.gamma}/ \
          -f -g {wildcards.gamma}"

rule gene_selection_for_correlation:
  input: "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/selectedGenes.txt"
  singularity: config["sif_file"]
  shell: "Rscript -e 'seurat.sc <- readRDS(\"{input}\");\
          selectedGenes <- rownames(seurat.sc[[\"RNA\"]]@counts)[rowSums(seurat.sc[[\"RNA\"]]@counts)> ncol(seurat.sc)*0.005];\
          write.table(selectedGenes,\"{output}\")'"

rule GA_computation_sc:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat_multimodal.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
  output: "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds"
  singularity: config["sif_file"]
  params: genome = config["pbmcMultiome"]["genome"]
  shell: "Rscript R/compute_activities.R -i {input.singlecells} \
          -o output/pbmcMultiome/singlecells_analysis/ \
          -w {threads} -g {params.genome} -u 2000 -d 100 -a ATAC -r -c -p {input.motifs}"

rule GA_computation_mc:
  input:
        metacells = "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mc.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
  output: "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds"
  singularity: config["sif_file"]
  params: genome = config["pbmcMultiome"]["genome"]
  shell: "Rscript R/compute_activities.R -i {input.metacells} \
          -o output/pbmcMultiome/{wildcards.inputMetacells}/g{wildcards.gamma}/ \
          -w {threads} -g {params.genome} -u 2000 -d 100 -a ATAC -r -c -p {input.motifs}"

rule mc_metrics_computation:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
        metacells = "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds"
  output: "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mcMetrics.rds"
  singularity: config["sif_file"]
  params:
        outdir = "output/pbmcMultiome/{inputMetacells}/g{gamma}/",
        python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python"
  shell: "Rscript R/compute_MC_metrics.R -s {input.singlecells} -o {params.outdir} -m {input.metacells} -n 1:50 -q 2:50 -e {params.python}"





# rule metacell_identification_pbmc_seacellsRNA_testRepro:
#   input:
#         singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.RNA.h5ad"
#   output: "output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/seacellMemberships.csv"
#   singularity: config["sif_file"]
#   shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
#           -o output/pbmcMultiome/test_seacellsRNA_rep/g{wildcards.gamma}/rep{wildcards.rep}/ \
#           -d 1:40 -r pca -g {wildcards.gamma}"
#
# rule data_aggregation_pbmc_seacellsRNA_testRepro:
#   input:
#         singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds",
#         memberships = "output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/seacellMemberships.csv"
#   output: "output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/seurat.multiome.mc.rds"
#   singularity: config["sif_file"]
#   params: workdir = wdir
#   shell: "export PATH={params.workdir}/config/bin/htslib-1.16/bin:$PATH;\
#           Rscript R/SCimplify10xMultiome_v5_CL.R -i {input.singlecells} -c {input.memberships} \
#           -o output/pbmcMultiome/test_seacellsRNA_rep/g{wildcards.gamma}/rep{wildcards.rep}/ \
#           -f -g {wildcards.gamma} -x SEACell-"
# rule mc_metrics_computation_testRepro:
#   input:
#         singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
#         metacells = "output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/seurat.multiome.mc.rds"
#   output: "output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/repro.seurat.multiome.mcMetrics.rds"
#   singularity: config["sif_file"]
#   params:
#         outdir = "output/pbmcMultiome/test_seacellsRNA_rep/g{gamma}/rep{rep}/",
#         python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python"
#   shell: "Rscript R/compute_MC_metrics.R -s {input.singlecells} -o {params.outdir} -m {input.metacells} -n 1:50 -q 2:50 -e {params.python}; mv {params.outdir}seurat.multiome.mcMetrics.rds {params.outdir}repro.seurat.multiome.mcMetrics.rds"

####################################################################################################################################
################################################ Test semi-supervised  #############################################################
####################################################################################################################################



rule Semi_sup_test_pbmc_multiome:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv"
  singularity: config["sif_file"]
  params: workdir = wdir
  # benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{graining}/seurat.multiome.mc.txt"
  shell: "Rscript R/SCimplify_test_semi_sup_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellMulti/testSemiSup/g{wildcards.graining}/ \
          -p 1:40 -q 2:40 -r SCT -a ATAC -v pca -w lsi -g {wildcards.graining} -s 2025 -l seurat_annotations"



####################################################################################################################################
################################################ Test supervised with noise  #######################################################
####################################################################################################################################



rule Shuffling_sup_test_pbmc_multiome:
  input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellMulti/testShufflingSup/g{graining}/results_test_shuffle_closest.csv"
  singularity: config["sif_file"]
  params: workdir = wdir
  # benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{graining}/seurat.multiome.mc.txt"
  shell: "Rscript R/SCimplify_test_shuffle_closest_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellMulti/testShufflingSup/g{wildcards.graining}/ \
          -p 1:40 -q 2:40 -r SCT -a ATAC -v pca -w lsi -g {wildcards.graining} -s 2025 -l seurat_annotations"


####################################################################################################################################
################################################ Test k_nn ################# #######################################################
####################################################################################################################################



rule k_nn_test_pbmc_multiome:
  input:
        script = "R/SCimplify_test_k_nn_CL.R",
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellMulti/testKNN/g{graining}/results_test_k_nn.csv"
  singularity: config["sif_file"]
  params: workdir = wdir
  # benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{graining}/seurat.multiome.mc.txt"
  shell: "Rscript R/SCimplify_test_k_nn_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellMulti/testKNN/g{wildcards.graining}/ \
          -p 1:40 -q 2:40 -r SCT -a ATAC -v pca -w lsi -g {wildcards.graining}"


rule k_nn_test_pbmc_multiome_sup:
  input:
        script = "R/SCimplify_test_k_nn_CL.R",
        singlecells = "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcMultiome/SuperCellMulti/testKNN_sup/g{graining}/results_test_k_nn.csv"
  singularity: config["sif_file"]
  params: workdir = wdir
  # benchmark :  "benchmark/pbmcMultiome/SuperCellMulti/g{graining}/seurat.multiome.mc.txt"
  shell: "Rscript R/SCimplify_test_k_nn_CL.R -i {input.singlecells} \
          -o output/pbmcMultiome/SuperCellMulti/testKNN_sup/g{wildcards.graining}/ \
          -p 1:40 -q 2:40 -r SCT -a ATAC -v pca -w lsi -g {wildcards.graining} -l seurat_annotations"




####################################################################################################################################
########################################## Benchmark and Correlations  #############################################################
####################################################################################################################################

rule atacRnaSinglecellAnalysis:
    input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
    output:
        "output/pbmcMultiome/singlecells_analysis/wilcox_rna.rds",
        "output/pbmcMultiome/singlecells_analysis/wilcox_motifs.rds",
        "output/pbmcMultiome/singlecells_analysis/multimodalMarkers_ttest.rds",
        "output/pbmcMultiome/singlecells_analysis/multimodalMarkers_wilcoxonAUC.rds",
        "output/pbmcMultiome/singlecells_analysis/CorrTables.rds"
    params:
        outdir = "output/pbmcMultiome/singlecells_analysis/",
        genome = config["pbmcMultiome"]["genome"],
        python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python"
    resources:
        mem_mb = 256000
    singularity: config["sif_file"]
    shell: "Rscript R/AtacRnaCorrAnalysis.R -i {input.singlecells} -o {params.outdir} -s {input.singlecells}"

rule atacRnaMetacellAnalysis:
    input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
        metacells = "output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
    output:
        "output/pbmcMultiome/{inputMetacells}/g{gamma}/multimodalMarkers_surveyweightedt.rds",
        "output/pbmcMultiome/{inputMetacells}/g{gamma}/CorrTables.rds"
    params:
        outdir = "output/pbmcMultiome/{inputMetacells}/g{gamma}/",
        genome = config["pbmcMultiome"]["genome"]
    resources:
        mem_mb = 100000
    singularity: config["sif_file"]
    shell: "Rscript R/AtacRnaCorrAnalysis.R -i {input.metacells} -o {params.outdir} -s {input.singlecells}"

rule MultimodalMarkers_mainCelltypes_sc:
    input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
        motifs =  "input/H13CORE_human_pfm.rds"
    output:
        "output/pbmcMultiome/singlecells_analysis/wilcox_rna_mainCellTypes.rds",
        "output/pbmcMultiome/singlecells_analysis/wilcox_motifs_mainCellTypes.rds",
        "output/pbmcMultiome/singlecells_analysis/multimodalMarkers_mainCellTypes_ttest.rds",
        "output/pbmcMultiome/singlecells_analysis/multimodalMarkers_mainCellTypes_wilcoxonAUC.rds"
    params:
        outdir = "output/pbmcMultiome/singlecells_analysis/",
        genome = config["pbmcMultiome"]["genome"]
    resources:
        mem_mb = 100000
    singularity: config["sif_file"]
    shell: "Rscript R/AtacRnaCorrAnalysis_mainCellTypes.R -i {input.singlecells} -o {params.outdir} -s {input.singlecells}"


# rule MultimodalMarkers_mainCelltypes:
#     input:
#         singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
#         metacells = "output/pbmcMultiome/SuperCellMulti/g{gamma}/seurat.multiome.activities.rds"
#     output:
#         "output/pbmcMultiome/SuperCellMulti/g{gamma}/multimodalMarkers_mainCellTypes_surveyweightedt.rds"
#     params:
#         outdir = "output/pbmcMultiome/SuperCellMulti/g{gamma}/",
#         genome = config["pbmcMultiome"]["genome"]
#     resources:
#         mem_mb = 100000
#     singularity: config["sif_file"]
#     shell: "Rscript R/AtacRnaCorrAnalysis_mainCellTypes.R -i {input.metacells} -o {params.outdir} -s {input.singlecells}"

rule pando_sc_pbmc:
    input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds"
    output:
        "output/pbmcMultiome/singlecells_analysis/grn_object_{pandoSC}.rds"
    params:
        outdir = "output/pbmcMultiome/singlecells_analysis/"
    benchmark: "benchmark/pbmcMultiome/singlecells_analysis/pando_{pandoSC}_pbmc.txt"
    resources:
        mem_mb = 64000
    singularity: config["sif_file"]
    shell: "Rscript R/pando_analysis.r -s {input.singlecells} -m {input.singlecells} -o {params.outdir} -w {threads} -x {wildcards.pandoSC}"

rule pando_mc_pbmc:
    input:
        singlecells = "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
        metacells = "output/pbmcMultiome/SuperCellMulti/g{gamma}/seurat.multiome.activities.rds",
        membership = "output/pbmcMultiome/SuperCellMulti/g{gamma}/SuperCellMemberships.csv"
    output:
        "output/pbmcMultiome/{inputMetacells}/g{gamma}/grn_object_{pandoMC}.rds"
    params:
        outdir = "output/pbmcMultiome/{inputMetacells}/g{gamma}/"
    benchmark: "benchmark/pbmcMultiome/{inputMetacells}/pando_mc_{pandoMC}_g{gamma}_pbmc.txt"
    resources:
        mem_mb = 64000
    singularity: config["sif_file"]
    shell: "Rscript R/pando_analysis.r -s {input.singlecells} -m {input.metacells} -b {input.membership} -o {params.outdir} -w {threads} -x {wildcards.pandoMC}"

rule gather_pando_pbmc:
    input:
        grn_res_mc = expand("output/pbmcMultiome/SuperCellMulti/g{gamma}/grn_object_{pandoMC}.rds", gamma = GAMMA, pandoMC = pando_method_mc),
        grn_res_sc = expand("output/pbmcMultiome/singlecells_analysis/grn_object_{pandoSC}.rds", pandoSC = pando_method_sc),
    output:
        "output/pbmcMultiome/SuperCellMulti/pando_metrics_summary_p_thresh0.1.RData"
    params:
        outdir = "output/pbmcMultiome/SuperCellMulti/"
    resources:
        mem_mb = 64000
    singularity: config["sif_file"]
    shell: "Rscript R/gather_pando_results.r -i 'output/pbmcMultiome/' -p 0.1 -o {params.outdir}"

rule gather_pando_pbmc2:
    input:
        grn_res_mc = expand("output/pbmcMultiome/SuperCellMulti/g{gamma}/grn_object_{pandoMC}.rds", gamma = GAMMA, pandoMC = pando_method_mc),
        grn_res_sc = expand("output/pbmcMultiome/singlecells_analysis/grn_object_{pandoSC}.rds", pandoSC = pando_method_sc),
    output:
        "output/pbmcMultiome/SuperCellMulti/pando_metrics_summary_p_thresh0.05.RData"
    params:
        outdir = "output/pbmcMultiome/SuperCellMulti/"
    resources:
        mem_mb = 64000
    singularity: config["sif_file"]
    shell: "Rscript R/gather_pando_results.r -i 'output/pbmcMultiome/' -p 0.05 -o {params.outdir}"

####################################################################################################################################
########################################## Figures #################################################################################
####################################################################################################################################

rule generate_figures_from_figure2:
  input:
    "figures/manuscript/final_figures_Rmd/Figure2_pbmcMultiome_bench_v2.Rmd",
    "output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds",
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.mcMetrics.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsMOFA","metaqMulti","metaqRNA","metaqATAC","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/seurat.multiome.activities.rds", gamma = GAMMA, inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsMOFA","metaqMulti","metaqRNA","metaqATAC","seacellsRNA","seacellsATAC","MetaCellRNA"]),
    expand("output/pbmcMultiome/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv", graining = ["20","75"]),
    expand("output/pbmcMultiome/SuperCellMulti/testShufflingSup/g{graining}/results_test_shuffle_closest.csv", graining = ["20","75"]),
    expand("output/pbmcMultiome/SuperCellMulti/testKNN/g{graining}/results_test_k_nn.csv", graining = ["20","75"]),
    expand("output/pbmcMultiome/SuperCellMulti/testKNN_sup/g{graining}/results_test_k_nn.csv", graining = ["20","75"])
  output: "figures/manuscript/final_figures_Rmd/Figure2_pbmcMultiome_bench_v2.html"
  singularity: config["sif_file"]
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"

# rule report_all_bench_corr:
#   input: "reports/pbmcMultiome/pbmcMultiome.Rmd",
#         "output/pbmcMultiome/singlecells_analysis/seurat.multiome.mc.activities.rds",
#         expand("output/pbmcMultiome/SuperCellMulti/testSemiSup/g{graining}/results_test_semisup.csv", graining = ["20","75"]),
#         expand("output/pbmcMultiome/{inputMetacells}/g{gamma}/corrTablePearson.csv",gamma = GAMMA,inputMetacells = ["SuperCellMulti","randomMetacells","SuperCellATAC","SuperCellRNA","seacellsRNA","seacellsATAC","MetaCellRNA"])
#   output: "reports/pbmcMultiome/pbmcMultiome.html"
#   container: "config/matk_multiomics.sif"
#   shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
