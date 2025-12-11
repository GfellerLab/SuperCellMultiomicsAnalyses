#############################################
##########  PBMC CITE-seq ATLAS   ###########
#############################################

rule download_pbmc_cite_seq_atlas:
  output:
    "input/GSM5008737_RNA_3P/matrix.mtx.gz",
    "input/GSM5008737_RNA_3P/features.tsv.gz",
    "input/GSM5008738_ADT_3P/matrix.mtx.gz",
    "input/GSM5008738_ADT_3P/features.tsv.gz"
    # "input/GSE164378_sc.meta.data_3P.csv"
  shell:
        """
        cd input;
        wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz" -O GSE164378_sc.meta.data_3P.csv.gz;
        gunzip GSE164378_sc.meta.data_3P.csv.gz;
        wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file" -O GSE164378_RAW.tar;
        mkdir -p GSE164378_RAW; tar -C GSE164378_RAW -xvf GSE164378_RAW.tar;
        cd GSE164378_RAW;
        mkdir GSM5008737_RNA_3P; mv GSM5008737_RNA_3P-* GSM5008737_RNA_3P/;
        mkdir GSM5008738_ADT_3P; mv GSM5008738_ADT_3P-* GSM5008738_ADT_3P/;
        cd GSM5008737_RNA_3P/;
        for file in GSM5008737_RNA_3P-*;
        do
          mv "$file" "${{file#GSM5008737_RNA_3P-}}";
        done;
        cd ../GSM5008738_ADT_3P/;
        for file in GSM5008738_ADT_3P-*;
        do
          mv "$file" "${{file#GSM5008738_ADT_3P-}}";
        done;
        cd ../../..;
        rm -f input/GSE164378_RAW/*.gz;
        mv input/GSE164378_RAW/* input/;
        rm -r input/GSE164378_RAW
        """
rule download_adt_table_pbmc_cite_seq_atlas:
  output: "input/pbmcCiteSeqAtlas/1-s2.0-S0092867421005833-mmc1.xlsx"
  shell: "cd input/pbmcCiteSeqAtlas/;wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867421005833-mmc1.xlsx"



rule load_pbmc_cite_seq_atlas:
  input: "input/GSM5008737_RNA_3P/matrix.mtx.gz",
          "input/GSM5008737_RNA_3P/features.tsv.gz",
          "input/GSM5008738_ADT_3P/matrix.mtx.gz",
          "input/GSM5008738_ADT_3P/features.tsv.gz"
          # "input/GSE164378_sc.meta.data_3P.csv",
  output: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  singularity: sif_file
  shell: "Rscript -e 'rnaCounts <- Seurat::Read10X(\"input/GSM5008737_RNA_3P/\");\
                      adtCounts <- Seurat::Read10X(\"input/GSM5008738_ADT_3P/\");\
                      meta <- read.csv(\"input/GSE164378_sc.meta.data_3P.csv\",row.names = 1);\
                      pbmc <- Seurat::CreateSeuratObject(counts = rnaCounts,meta.data = meta);\
                      pbmc[[\"ADT\"]] <- Seurat::CreateAssayObject(counts = adtCounts);\
                      pbmc <- pbmc[,pbmc$celltype.l2 != \"Doublet\"];\
                      pbmc$celltype.l1.5 <- pbmc$celltype.l1;\
                      pbmc$celltype.l1.5[pbmc$celltype.l1 == \"other\"] <- pbmc$celltype.l2[pbmc$celltype.l1 == \"other\"];\
                      pbmc$celltype.l1.5[pbmc$celltype.l1 == \"other T\"] <- pbmc$celltype.l2[pbmc$celltype.l1 == \"other T\"];\
                      saveRDS(pbmc,\"{output}\")'"

rule one_sample_pbmc_cite_atlas_wnn_analyzis:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds"
  params: smpName = "{pbmcCiteSmp}",
          outdir = "output/pbmcCiteSeqAtlas/{params.pbmcCiteSmp}/singlecells_analysis",
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnAnalysisPBMC_CiteAtlasCL.R -i {input} \
         -s {wildcards.pbmcCiteSmp} \
         -r SCT \
         -d -y {params.python} \
         -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/singlecells_analysis \
         -p 1:40 -q 1:50"

rule singlecells_rna_pca_analysis_one_sample_pbmc_cite_atlas:
  input:  sobj =  "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.h5ad",
  singularity: sif_file
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/RNA_pca_CiteSeqCL.R -i {input.sobj} \
          -r SCT -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/singlecells_analysis"

rule singlecells_adt_pca_analysis_one_sample_pbmc_cite_atlas:
  input: sobj =  "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.ADT.h5ad",
  singularity: sif_file
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.ADT.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/ADT_pca_CiteSeqCL.R -i {input.sobj} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/singlecells_analysis"

rule preprocessing_supercell_one_sample_pbmc_cite_atlas:
  input: sobj =  "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat_multimodal.rds",
  singularity: sif_file
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat_multimodal.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/preprocessing_seurat_for_SuperCell_CiteSeqCL.R -i {input.sobj} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/singlecells_analysis "

rule metacell_identification_sample_pbmc_cite_atlas_seacellsRNA:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/seacellMemberships.csv"
  singularity: sif_file
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/seacellMemberships.txt"
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsRNA/g{wildcards.gamma}/\
          -d 1:40 -r pca -g {wildcards.gamma}"


rule metacell_identification_sample_pbmc_cite_atlas_MetaCellRNA:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/MetaCellMemberships.txt"
  singularity: sif_file
  shell: "python3 python/MATK_MetaCell2CL.py -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/MetaCellRNA/g{wildcards.gamma}/\
          -g {wildcards.gamma}"



rule metacell_identification_sample_pbmc_cite_atlas_seacellsADT:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.ADT.h5ad",
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/seacellMemberships.csv"
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/seacellMemberships.txt"
  singularity: sif_file
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsADT/g{wildcards.gamma}/ \
          -d 1:50 -r apca -g {wildcards.gamma}"



rule metacell_identification_pbmc_cite_atlas:
  input:
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat_multimodal.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/SuperCellHierarchy.rds"
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/SuperCellHierarchy.txt"
  singularity: sif_file
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:40 -q 1:50 -v 1:40 -w 1:50 -r SCT -e TRUE -k 30 -g {wildcards.gamma}\
          -t {params.python} -s SuperCellHierarchy"

rule data_aggregation_pbmc_cite_atlas:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        hierarchy = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/SuperCellHierarchy.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/seurat.cite.mc.rds",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -m {input.hierarchy} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellMulti/g{wildcards.gamma}/ \
          -g {wildcards.gamma} \
          -y SuperCell_Multi \
          -z {input.gene_protein}\
          -t {params.python}\
          -s seurat"

rule metacell_identification_pbmc_cite_atlas_ADT:
  input:
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.ADT.h5ad"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/SuperCellHierarchy.rds"
  singularity: sif_file
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/SuperCellHierarchy.txt"
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellADT/g{wildcards.gamma}/ \
          -q 1:50 -v 1:40 -w 1:50 -r SCT -e TRUE -k 30 -g {wildcards.gamma}\
          -t {params.python} -s SuperCellHierarchy"

rule data_aggregation_pbmc_cite_atlas_ADT:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        hierarchy = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/SuperCellHierarchy.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/seurat.cite.mc.rds",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/SuperCellHierarchy.txt"
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -m {input.hierarchy} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellADT/g{wildcards.gamma}/ \
          -g {wildcards.gamma} \
          -y SuperCell_ADT \
          -z {input.gene_protein}\
          -t {params.python}\
          -s seurat"

rule metacell_identification_pbmc_cite_atlas_RNA:
  input:
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/SuperCellHierarchy.rds"
  singularity: sif_file
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  benchmark: "benchmark/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/SuperCellHierarchy.txt"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellRNA/g{wildcards.gamma}/ \
          -p 1:40 -v 1:40 -w 1:50 -r SCT -e TRUE -k 30 -g {wildcards.gamma}\
          -t {params.python} -s SuperCellHierarchy"

rule data_aggregation_pbmc_cite_atlas_RNA:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        hierarchy = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/SuperCellHierarchy.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/seurat.cite.mc.rds",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  params: workdir = wdir,
          python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -m {input.hierarchy} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} \
          -y SuperCell_RNA \
          -z {input.gene_protein}\
          -t {params.python}\
          -s seurat"






rule random_metacell_pbmc_cite_atlas:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/randomMetacells/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/randomMetacells/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/randomMetacells/g{wildcards.gamma}/ \
          -d TRUE -v 1:40 -w 1:50 -e TRUE -g {wildcards.gamma}\
          -y randomMetacells -z {input.gene_protein}\
          -t {params.python}"


rule data_aggregation_pbmc_cite_atlas_seacellsRNA:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        memberships = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -v 1:40 -w 1:50 \
          -y SEACells_RNA -z {input.gene_protein}\
          -t {params.python}"

rule data_aggregation_pbmc_cite_MetaCellRNA:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        memberships = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -v 1:40 -w 1:50 \
          -y MetaCell_RNA -z {input.gene_protein}\
          -t {params.python}"

rule data_aggregation_pbmc_cite_atlas_seacellsADT:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        memberships = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/metaData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/corrTablePearson.csv"
  singularity: sif_file
  params: python = "/opt/conda/envs/MetacellAnalysisToolkit/bin/python3.9"
  shell: "Rscript R/SCimplifyCiteSeq_v5_CL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsADT/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -v 1:40 -w 1:50 \
          -y SEACells_ADT -z {input.gene_protein}\
          -t {params.python}"


### integration single cell CITE-atlas all samples ###

# rule cite_atlas_integrated_wnn_analyzis:
#   input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells/seuratCombinedWNN.rds"
#   singularity: sif_file
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells/bench.txt"
#   shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasCL.R -i {input} \
#          -r SCT \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells/ \
#          -p 1:40 -q 1:50"

rule cite_atlas_singlecells_SCT_seuratRPCA_SCT:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/bench.txt"
  shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasScLevelCL.R -i {input} \
         -r SCT \
         -o output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/ \
         -p 1:40 -q 1:50"

rule bench_cite_atlas_singlecells_SCT_seuratRPCA_SCT:
  input: "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/bench_res.rds"
  singularity: sif_file
  shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_single_cells_CL.R -i {input} \
         -r SCT \
         -o output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/ \
         -p 1:40 -q 1:50"

# rule cite_atlas_integrated_wnn_analyzis_RNA_only:
#   input: single_cells = "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds",
#          install =  "config/installed_sup_tools"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_unsup/seuratCombinedWNN.rds"
#   conda: SupervisedTools
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_unsup/seuratCombinedWNN.txt"
#   shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasMcLevelCL_RNA_unsup.R -i {input.single_cells} \
#          -g 20 \
#          -m SCT \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_unsup/ \
#          -p 1:40 -q 1:50"
#
# rule cite_atlas_integrated_wnn_analyzis_RNA_scATOMIC:
#   input: single_cells = "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds",
#          install =  "config/installed_sup_tools"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scATOMIC/seuratCombinedWNN.rds"
#   conda: SupervisedTools
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scATOMIC/seuratCombinedWNN.txt"
#   shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasSupMcLevelCL_RNA.R -i {input.single_cells} \
#          -g 20 \
#          -c scATOMIC \
#          -m SCT \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scATOMIC/ \
#          -p 1:40 -q 1:50"
#
# rule cite_atlas_integrated_wnn_analyzis_RNA_scGate:
#   input: single_cells = "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds",
#          install =  "config/installed_sup_tools"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scGate/seuratCombinedWNN.rds"
#   conda: SupervisedTools
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scGate/seuratCombinedWNN.txt"
#   shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasSupMcLevelCL_RNA.R -i {input.single_cells} \
#          -g 20 \
#          -c scGate \
#          -m SCT \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/RNA_scGate/ \
#          -p 1:40 -q 1:50"
#
#
# rule cite_monoCD14_atlas_integrated_wnn_analyzis:
#   input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells_MonoCD14/seuratCombinedWNN.rds"
#   singularity: sif_file
#   benchmark:
#     "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells_MonoCD14/bench.txt"
#   shell: "Rscript  snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisMonoCD14_CiteAtlasCL.R -i {input} \
#          -r SCT \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells_MonoCD14/ \
#          -p 1:40 -q 1:50"
#
# rule cite_atlas_integrated_wnn_analyzis_MC_logNorm:
#   input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g20/fullLogNorm/seuratCombinedWNN.rds"
#   singularity: sif_file
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g20/bench.txt"
#   shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasMcLevelCL.R -i {input} \
#          -g 20 \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g20/fullLogNorm \
#          -p 1:40 -q 1:50"
#
rule cite_atlas_unsupMetacells_SCT_seuratRPCA_SCT:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasMcLevelCL.R -i {input} \
         -g 20 \
         -r SCT \
         -m SCT \
         -o output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/ \
         -p 1:40 -q 1:50"

rule bench_cite_atlas_unsupMetacells_SCT_seuratRPCA_SCT:
  input: mc = "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc} \
         -r SCT \
         -o output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/ \
         -p 1:40 -q 1:50"

# rule cite_atlas_unsupMetacells_SCT_unsupSTACAS_lognorm:
#   input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
#   output: "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/seuratCombinedWNN.rds"
#   singularity: sif_file
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/bench.txt"
#   shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasMcLevelCL.R -i {input} \
#          -g 20 \
#          -r SCT \
#          -m SCT \
#          -o output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/ \
#          -p 1:40 -q 1:50"
#
# rule cite_atlas_unsupMetacells_SCT_unsupSTACAS_lognorm:
#   input: mc = "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/seuratCombinedWNN.rds",
#          sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
#   output: "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/bench_res.rds"
#   singularity: sif_file
#   shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
#          -j {input.sc} \
#          -r SCT \
#          -o output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/ \
#          -p 1:40 -q 1:50"
#
# rule cite_atlas_integrated_wnn_analyzis_Sup_MC_SCT_STACAS_logNorm:
#   input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
#   output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/seuratCombinedWNN.rds"
#   singularity: sif_file
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm/bench.txt"
#   shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasSupMcLevelCL.R -i {input} \
#          -g 20 \
#          -m SCT \
#          -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_sup_metacells_supSTACAS/g20/metacell_SCT_STACAS_logNorm \
#          -p 1:40 -q 1:50"


rule cite_atlas_100percent_semisupMetacells_SCT_supStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -e 1 \
         -o output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule bench_cite_atlas_100percent_semisupMetacells_SCT_supStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc} \
         -o output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule cite_atlas_50percent_semisupMetacells_SCT_supStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -e 0.50 \
         -o output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule bench_cite_atlas_50percent_semisupMetacells_SCT_supStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc} \
         -o output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule cite_atlas_20percent_semisupMetacells_SCT_supStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -e 0.20 \
         -o output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule bench_cite_atlas_20percent_semisupMetacells_SCT_supStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc}\
         -o output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule cite_atlas_5percent_semisupMetacells_SCT_supStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -e 0.05 \
         -o output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule bench_cite_atlas_5percent_semisupMetacells_SCT_supStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc}  \
         -o output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule cite_atlas_supMetacells_SCT_supStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -o output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"



rule bench_cite_atlas_supMetacells_SCT_supStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc} \
         -o output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"


rule cite_atlas_supMetacells_SCT_unsupStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -u \
         -o output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"



rule bench_cite_atlas_supMetacells_SCT_unsupStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc} \
         -o output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"


rule cite_atlas_unsupMetacells_SCT_supStacas_lognorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  singularity: sif_file
  benchmark:
    "benchmark/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/bench.txt"
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
         -g 20 \
         -m SCT \
         -v \
         -o output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"



rule bench_cite_atlas_unsupMetacells_SCT_supStacas_lognorm:
  input: mc = "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds",
         sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
  output: "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds"
  singularity: sif_file
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
         -j {input.sc} \
         -o output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/ \
         -p 1:40 -q 1:50"

rule run_edgeR_mono:
  input: mc = "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  output:
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_all.txt",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_CD14.txt",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/edgeR_res_pairwise.txt"
  singularity: sif_file_edger
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/edgeR_diff_analysis.R -c {input.mc} \
         -o output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/"

rule run_edgeR_mono_ADT:
  input: mc = "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/seuratCombinedWNN.rds"
  output:
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ADT_diff/edgeR_res_all.txt",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ADT_diff/edgeR_res_CD14.txt",
    "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ADT_diff/edgeR_res_pairwise.txt"
  singularity: sif_file_edger
  shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/edgeR_diff_analysis.R -c {input.mc} \
         -o output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/ADT_diff/ -a ADT"

# rule cite_atlas_supMetacells_SCT_supStacas_lognorm_mean:
#   input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
#   output: "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/seuratCombinedWNN.rds"
#   singularity: sif_file
#   benchmark:
#     "benchmark/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/bench.txt"
#   shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasSupMcLevelCL.R -i {input} \
#          -g 20 \
#          -d logMean \
#          -m SCT \
#          -o output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/ \
#          -p 1:40 -q 1:50"
#
# rule bench_cite_atlas_supMetacells_SCT_supStacas_lognorm_mean:
#   input: mc = "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/seuratCombinedWNN.rds",
#          sc = "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/seuratCombinedWNN.rds"
#   output: "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/bench_res.rds"
#   singularity: sif_file
#   shell: "Rscript snakemakeWorkflows/pbmcCiteSeqAtlas/R/CiteAtlasBenchmark_CL.R -i {input.mc} \
#          -j {input.sc} \
#          -o output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/ \
#          -p 1:40 -q 1:50"

# rule report_pbmc_cite_atlas_samples:
#   input: "reports/pbmcCiteSeqAtlas/Correlation_MC_metrics_pbmc_cite_atlas_analysis.Rmd",
#         "input/pbmcCiteSeqAtlas/1-s2.0-S0092867421005833-mmc1.xlsx",
#          expand("output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/{method}/g{gamma}/metaData.csv",pbmcCiteSmp= pbmcCiteSamples, gamma = GAMMA,method = CiteAtlasMethod),
#          expand("output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/{method}/g{gamma}/corrTablePearson.csv",pbmcCiteSmp= pbmcCiteSamples, gamma = GAMMA,method = CiteAtlasMethod)
#   output: "reports/pbmcCiteSeqAtlas/Correlation_MC_metrics_pbmc_cite_atlas_analysis.html"
#   singularity: sif_file
#   shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
#
# rule report_pbmc_cite_atlas_integration:
#   input: "reports/pbmcCiteSeqAtlas/Benchmark_pbmc_cite_atlas_integration_metacell_analysis.Rmd",
#          "output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/bench_res.rds",
#          "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
#          "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
#           # "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean/g20/bench_res.rds",
#          "output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
#          "output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
#          "output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
#           "output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
#          "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/bench_res.rds"
#          # "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_seuratRPCA_SCT/g20/bench_res.rds"
#   output: "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis.html"
#   singularity: sif_file
#   shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"


rule generate_figures_from_pbmc_cite_atlas:
  input: "figures/manuscript/final_figures_Rmd/Figure4_bench_pbmc_cite_atlas_integration_v2.Rmd",
         "output/pbmcCiteSeqAtlas/supMetacells_SCT_unsupStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/unsupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/5percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/20percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/50percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/100percent_semisupMetacells_SCT_supStacas_lognorm/g20/bench_res.rds",
         "output/pbmcCiteSeqAtlas/singlecells_SCT_seuratRPCA_SCT/bench_res.rds"
  output: "figures/manuscript/final_figures_Rmd/Figure4_bench_pbmc_cite_atlas_integration_v2.html"
  singularity: sif_file
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
