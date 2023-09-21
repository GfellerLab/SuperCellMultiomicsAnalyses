#############################################
##########  PBMC CITE-seq ATLAS   ###########
#############################################
          
rule download_pbmc_cite_seq_atlas:
  output: "input/GSE164378_sc.meta.data_3P.csv",
          "input/GSM5008737_RNA_3P/matrix.mtx.gz",
          "input/GSM5008737_RNA_3P/features.tsv.gz",
          "input/GSM5008737_ADT_3P/matrix.mtx.gz",
          "input/GSM5008737_ADT_3P/features.tsv.gz"
  shell: "cd input;\
          wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz;\
          gunzip GSE164378_sc.meta.data_3P.csv.gz;\
          wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file;\
          tar -C GSE164378_RAW -xvf GSE164378_RAW.tar;cd GSE164378_RAW;\
          mkdir GSM5008737_RNA_3P; mv GSM5008738_RNA_3P-*;\
          mkdir GSM5008738_ADT_3P; mv GSM5008738_ADT_3P-*;\
          cd ../input/GSM5008737_RNA_3P/;\
          for file in GSM5008737_RNA_3P-*;\
          do;\
          mv \"$file\" \"${file#GSM5008737_RNA_3P-}\";\
          done;\
          cd ../input/GSM5008737_ADT_3P/;\
          for file in GSM5008737_ADT_3P-*;\
          do;\
          mv \"$file\" \"${file#GSM5008737_ADT_3P-}\";\
          done;\
          cd ..;rm -f *.gz"
          
rule download_adt_table_pbmc_cite_seq_atlas:
  output: "input/pbmcCiteSeqAtlas/1-s2.0-S0092867421005833-mmc1.xlsx"
  shell: "cd input/pbmcCiteSeqAtlas/;wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867421005833-mmc1.xlsx"

          
          
rule load_pbmc_cite_seq_atlas:
  input: "input/GSE164378_sc.meta.data_3P.csv",
          "input/GSM5008737_RNA_3P/matrix.mtx.gz",
          "input/GSM5008737_RNA_3P/features.tsv.gz",
          "input/GSM5008738_ADT_3P/matrix.mtx.gz",
          "input/GSM5008738_ADT_3P/features.tsv.gz",
          "config/SuperCellMultiomics/installedSuperCellMultiomics"
  output: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'rnaCounts <- Seurat::Read10X(\"input/GSM5008737_RNA_3P/\");\
                      adtCounts <- Seurat::Read10X(\"input/GSM5008738_ADT_3P/\");\
                      meta <- read.csv(\"input/GSE164378_sc.meta.data_3P.csv\",row.names = 1);\
                      pbmc <- Seurat::CreateSeuratObject(counts = rnaCounts,meta.data = meta);\
                      pbmc[[\"ADT\"]] <- Seurat::CreateAssayObject(counts = adtCounts);\
                      pbmc <- pbmc[,pbmc$celltype.l2 != \"Doublet\"];\
                      saveRDS(pbmc,\"{output}\")'"

rule one_sample_pbmc_cite_atlas_wnn_analyzis:
  input: "output/pbmcCiteSeqAtlas/pbmc_hiv_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.SCT.h5ad",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.ADT.h5ad",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.h5ad"
  params: smpName = "{pbmcCiteSmp}",
          outdir = "output/pbmcCiteSeqAtlas/{params.pbmcCiteSmp}/singlecells_analysis"
  conda: SuperCellMultiomics
  shell: "Rscript R/wnnAnalysisPBMC_CiteAtlasCL.R -i {input} \
         -s {wildcards.pbmcCiteSmp} \
         -r SCT \
         -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/singlecells_analysis \
         -p 1:40 -q 1:50"
         
rule metacell_identification_sample_pbmc_cite_atlas_seacellsRNA:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.SCT.h5ad"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/seacellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsRNA/g{wildcards.gamma}/\
          -d 1:40 -r pca -g {wildcards.gamma}"
          
          
rule metacell_identification_sample_pbmc_cite_atlas_MetaCellRNA:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.RNA.h5ad"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  shell: "python3 python/MetaCell2CL.py -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/MetaCellRNA/g{wildcards.gamma}/\
          -g {wildcards.gamma}"

          

rule metacell_identification_sample_pbmc_cite_atlas_seacellsADT:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seurat.ADT.h5ad",
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/seacellMemberships.csv"
  conda: SuperCellMultiomics_pyenv
  shell: "python3 python/SEACellsCL.py -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsADT/g{wildcards.gamma}/ \
          -d 1:50 -r apca -g {wildcards.gamma}"
          
        
        
          
          
rule metacell_identification_pbmc_cite_atlas:
  input: 
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellMulti/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellMulti/g{wildcards.gamma}/ \
          -p 1:40 -q 1:50 -v 1:40 -w 1:50 -r SCT -e TRUE -k 30 -g {wildcards.gamma}\
          -y SuperCell_Multi -z {input.gene_protein} \
          -t {params.python}"
          
rule metacell_identification_pbmc_cite_atlas_RNA:
  input: 
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellRNA/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellRNA/g{wildcards.gamma}/ \
          -p 1:40 -v 1:40 -w 1:50 -r SCT -e TRUE -k 30 -g {wildcards.gamma} \
          -y SuperCell_RNA -z {input.gene_protein} \
          -t {params.python}"
          

          
rule metacell_identification_pbmc_cite_atlas_ADT:
  input:
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/SuperCellADT/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/SuperCellADT/g{wildcards.gamma}/ \
          -q 1:50 -v 1:40 -w 1:50 -e TRUE -k 30 -g {wildcards.gamma}\
          -y SuperCell_ADT -z {input.gene_protein}\
          -t {params.python}"
          
rule random_metacell_pbmc_cite_atlas:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/randomMetacells/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/randomMetacells/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/randomMetacells/g{wildcards.gamma}/ \
          -d TRUE -v 1:40 -w 1:50 -e TRUE -g {wildcards.gamma}\
          -y randomMetacells -z {input.gene_protein}\
          -t {params.python}"
          
          
rule data_aggregation_pbmc_cite_atlas_seacellsRNA:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        memberships = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsRNA/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -v 1:40 -w 1:50 \
          -y SEACells_RNA -z {input.gene_protein}\
          -t {params.python}"

rule data_aggregation_pbmc_cite_MetaCellRNA:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        memberships = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/MetaCellMemberships.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/MetaCellRNA/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/MetaCellRNA/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -v 1:40 -w 1:50 \
          -y MetaCell_RNA -z {input.gene_protein}\
          -t {params.python}"
          
rule data_aggregation_pbmc_cite_atlas_seacellsADT:
  input: 
        singlecells = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/singlecells_analysis/seuratWNN.rds",
        gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
        memberships = "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/seacellMemberships.csv"
  output: "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/crData.csv",
          "output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/seacellsADT/g{gamma}/metaData.csv"
  conda: SuperCellMultiomics
  params: python = Seacells + "/bin/python"
  shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} -c {input.memberships} \
          -o output/pbmcCiteSeqAtlas/{wildcards.pbmcCiteSmp}/seacellsADT/g{wildcards.gamma}/ \
          -g {wildcards.gamma} -v 1:40 -w 1:50 \
          -y SEACells_ADT -z {input.gene_protein}\
          -t {params.python}"
          
          
          
### integration single cell CITE-atlas ###

rule cite_atlas_integrated_wnn_analyzis:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells/seuratCombinedWNN.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/cite_atlas_integrated_wnn_analyzis.txt"
  shell: "Rscript  pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasCL.R -i {input} \
         -r SCT \
         -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells/ \
         -p 1:40 -q 1:50"
         
rule cite_atlas_integrated_wnn_analyzis_MC_logNorm:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g{gamma_atlas}/fullLogNorm/seuratCombinedWNN.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/cite_atlas_integrated_wnn_analyzis_MC/g{gamma_atlas}/fullLogNorm/cite_atlas_integrated_wnn_analyzis_MC.txt"
  shell: "Rscript pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasMcLevelCL.R -i {input} \
         -g {wildcards.gamma_atlas} \
         -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g{wildcards.gamma_atlas}/fullLogNorm \
         -p 1:40 -q 1:50"
         
rule cite_atlas_integrated_wnn_analyzis_MC_SCT:
  input: "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
  output: "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g{gamma_atlas}/fullSCT/seuratCombinedWNN.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/cite_atlas_integrated_wnn_analyzis_MC/g{gamma_atlas}/fullSCT/cite_atlas_integrated_wnn_analyzis_MC.txt"
  shell: "Rscript pbmcCiteSeqAtlas/R/wnnIntegratedAnalysisPBMC_CiteAtlasMcLevelCL.R -i {input} \
         -g {wildcards.gamma_atlas} \
         -r SCT \
         -m SCT \
         -o output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_metacells/g{wildcards.gamma_atlas}/fullSCT \
         -p 1:40 -q 1:50"
         
# rule SuperCellMulti_on_integrated_data:
#   input: 
#         gene_protein = "input/pbmcCiteSeqAtlas/gene_protein.csv",
#         singlecells = "output/pbmcCiteSeqAtlas/pbmc_HIV_integrated/seuratCombinedWNN.rds"
#   output: "output/pbmcCiteSeqAtlas/pbmc_HIV_integrated/SuperCellMulti/g50/crData.csv",
#           "output/pbmcCiteSeqAtlas/pbmc_HIV_integrated/SuperCellMulti/g50/metaData.csv"
#   conda: SuperCellMultiomics
#   params: python = Seacells + "/bin/python"
#   shell: "Rscript R/SCimplifyCiteSeqCL.R -i {input.singlecells} \
#           -o output/pbmcCiteSeqAtlas/pbmc_HIV_integrated/SuperCellMulti/g50/ \
#           -p 1:40 -q 1:50 -v 1:40 -w 1:50 -r SCT -e TRUE -k 30 -g 50\
#           -y SuperCell_Multi -z {input.gene_protein} \
#           -t {params.python}"

          
rule report_pbmc_cite_atlas_samples:
  input: "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_analysis.Rmd",
        "input/pbmcCiteSeqAtlas/1-s2.0-S0092867421005833-mmc1.xlsx",
         expand("output/pbmcCiteSeqAtlas/{pbmcCiteSmp}/{method}/g{gamma}/metaData.csv",pbmcCiteSmp= pbmcCiteSamples, gamma = GAMMA,method = CiteAtlasMethod),
  output: "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_analysis.html"
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
  
rule report_pbmc_cite_atlas_integration:
  input: "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_integration_metacell_analysis",
        "input/pbmcCiteSeqAtlas/1-s2.0-S0092867421005833-mmc1.xlsx",
         "output/correlationAnalyzis/CITEseq/pbmc_HIV_integrated_MC/g20/fullSCT//seuratCombinedWNN.rds"
  output: "reports/pbmcCiteSeqAtlas/pbmc_cite_atlas_analysis.html"
  shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
