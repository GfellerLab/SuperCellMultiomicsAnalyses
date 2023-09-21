##########################################################
######### Multiomic single cell velocities ###############
##########################################################

rule download_processed_data_multiomic_velocity:
  output: rna = "input/multiomicVelocity/{velocitySmp}/adata_postpro_rna.h5ad",
          atac = "input/multiomicVelocity/{velocitySmp}/adata_postpro_atac.h5ad"
  params: 
          rna = lambda w: config[w.velocitySmp]['processedRNA'],
          atac = lambda w: config[w.velocitySmp]['processedATAC']
  shell: "wget -O {output.rna} {params.rna} && wget -O {output.atac} {params.atac}"

rule download_raw_data_velocity_humanHSPC:
  input: "input/multiomicVelocity/humanHSPC/adata_postpro_rna.h5ad"
  output: rna = "input/multiomicVelocity/humanHSPC/adata_raw_rna.h5ad",
          atac = "input/multiomicVelocity/humanHSPC/adata_raw_atac.h5ad"
  params: mtx = config["humanHSPC"]['rawMatrices'],
          features = config["humanHSPC"]['rawFeatures'],
          cells = config["humanHSPC"]['rawCells'],
          clusterLabel = config["humanHSPC"]['label']
  run: 
    import wget
    try:
        os.mkdir("input/multiomicVelocity/humanHSPC/filtered_feature_bc_matrix/")
    except FileExistsError:
        pass
    wget.download(params.features, "input/multiomicVelocity/humanHSPC/filtered_feature_bc_matrix/features.tsv.gz")
    wget.download(params.mtx, "input/multiomicVelocity/humanHSPC/filtered_feature_bc_matrix/matrix.mtx.gz")
    wget.download(params.cells, "input/multiomicVelocity/humanHSPC/filtered_feature_bc_matrix/barcodes.tsv.gz")
    adata = sc.read_10x_mtx('input/multiomicVelocity/humanHSPC/filtered_feature_bc_matrix/', var_names='gene_symbols', gex_only=False)
    adata_atac = adata[:,adata.var['feature_types'] == "Peaks"]
    adata_rna = adata[:,adata.var['feature_types'] == "Gene Expression"]
    adata = sc.read_h5ad(input[0])
    adata_rna = adata_rna[adata.obs_names,]
    adata_rna.obs[params.clusterLabel] = adata.obs[params.clusterLabel]
    adata_atac = adata_atac[adata.obs_names]
    adata_rna.write(output.rna)
    adata_atac.write(output.atac) 
  

rule download_raw_data_velocity_humanBrain:
  input: "input/multiomicVelocity/humanBrain/adata_postpro_rna.h5ad"
  output: rna = "input/multiomicVelocity/humanBrain/adata_raw_rna.h5ad",
          atac = "input/multiomicVelocity/humanBrain/adata_raw_atac.h5ad"
  params: rawATAC = config["humanBrain"]['rawSceATAC'],
          rawRNA = config["humanBrain"]['rawSceRNA'],
          clusterLabel = config["humanBrain"]['label']
  conda: SuperCellMultiomics
  shell: "wget -O input/multiomicVelocity/humanBrain/Multiome_RNA_SCE.RDS {params.rawRNA};\
          wget -O input/multiomicVelocity/humanBrain/Multiome_ATAC_SCE.RDS {params.rawATAC};\
          Rscript -e 'library(SingleCellExperiment);sce <- readRDS(\"input/multiomicVelocity/humanBrain/Multiome_RNA_SCE.RDS\");\
          metaFiltered <- anndata::read_h5ad(\"{input[0]}\")$obs;\
          filteredCells <- rownames(metaFiltered);\
          cluster <- metaFiltered[,\"{params.clusterLabel}\"];\
          adata <- anndata::AnnData(X = t(assay(sce)),obs = data.frame(colData(sce)),var = data.frame(rowData(sce)));\
          adata <- adata[filteredCells,];\
          adata$obs[,\"{params.clusterLabel}\"] <- metaFiltered[,\"{params.clusterLabel}\"];\
          anndata::write_h5ad(adata,\"{output.rna}\");\
          sce <- readRDS(\"input/multiomicVelocity/humanBrain/Multiome_ATAC_SCE.RDS\");\
          adata <- anndata::AnnData(X = t(assay(sce)),obs = data.frame(colData(sce)),var = data.frame(rowData(sce)));\
          adata <- adata[filteredCells,];\
          adata$obs[,\"{params.clusterLabel}\"] <- metaFiltered[,\"{params.clusterLabel}\"];\
          anndata::write_h5ad(adata,\"{output.atac}\")'"
          
        
rule SuperCell_for_multivelo:
  input: rna = "input/multiomicVelocity/{velocitySmp}/adata_postpro_rna.h5ad",
         atac = "input/multiomicVelocity/{velocitySmp}/adata_postpro_atac.h5ad"
  output: "output/multiomicVelocity/{velocitySmp}/SuperCellPostPro/adata_mc_rna.h5ad",
          "output/multiomicVelocity/{velocitySmp}/SuperCellPostPro/adata_mc_atac.h5ad",
          "output/multiomicVelocity/{velocitySmp}/SuperCellPostPro/nn_idx.txt",
          "output/multiomicVelocity/{velocitySmp}/SuperCellPostPro/nn_dist.txt",
          "output/multiomicVelocity/{velocitySmp}/SuperCellPostPro/nn_cells.txt"
  params: clusterLabel =  lambda w: config[w.velocitySmp]['label'],
          cellCycleCr = lambda w: config[w.velocitySmp]['cellCycleOpt']
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/multiomicVelocity/{velocitySmp}/SuperCellPostPro_benchmark.txt"
  shell: "Rscript snakemakeWorkflows/multiomicVelocity/R/SuperCellPostProMultiveloCL.R \
  -r {input.rna} -a {input.atac} \
  -i RNA -g 10 -l {params.clusterLabel} \
  -o output/multiomicVelocity/{wildcards.velocitySmp}/SuperCellPostPro/ {params.cellCycleCr} -u"

rule Subsampling_for_multivelo:
  input: rnaPostPro = "input/multiomicVelocity/{velocitySmp}/adata_postpro_rna.h5ad",
         atacPostPro = "input/multiomicVelocity/{velocitySmp}/adata_postpro_atac.h5ad",
         rnaRaw = "input/multiomicVelocity/{velocitySmp}/adata_raw_rna.h5ad",
         atacRaw = "input/multiomicVelocity/{velocitySmp}/adata_raw_atac.h5ad"
  output: "output/multiomicVelocity/{velocitySmp}/Subsampling/adata_mc_rna.h5ad",
          "output/multiomicVelocity/{velocitySmp}/Subsampling/adata_mc_atac.h5ad",
          "output/multiomicVelocity/{velocitySmp}/Subsampling/nn_idx.txt",
          "output/multiomicVelocity/{velocitySmp}/Subsampling/nn_dist.txt",
          "output/multiomicVelocity/{velocitySmp}/Subsampling/nn_cells.txt"
  params: clusterLabel =  lambda w: config[w.velocitySmp]['label'],
          cellCycleCr = lambda w: config[w.velocitySmp]['cellCycleOpt']
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/multiomicVelocity/{velocitySmp}/Subsampling_benchmark.txt"
  shell: "Rscript snakemakeWorkflows/multiomicVelocity/R/SubsampleForMultiveloCL.R\
  -r {input.rnaPostPro} -a {input.atacPostPro}\
  -t {input.rnaRaw} -b {input.atacRaw}\
  -g 10 -l {params.clusterLabel} -o output/multiomicVelocity/{wildcards.velocitySmp}/Subsampling {params.cellCycleCr} -u"


rule accelerated_multivelo:
  input: rna =  "output/multiomicVelocity/{velocitySmp}/{runMv}/adata_mc_rna.h5ad",
         atac = "output/multiomicVelocity/{velocitySmp}/{runMv}/adata_mc_atac.h5ad",
         idx = "output/multiomicVelocity/{velocitySmp}/{runMv}/nn_idx.txt",
         dist = "output/multiomicVelocity/{velocitySmp}/{runMv}/nn_dist.txt",
         cells = "output/multiomicVelocity/{velocitySmp}/{runMv}/nn_cells.txt"
  output: "output/multiomicVelocity/{velocitySmp}/{runMv}/multivelo_result.h5ad"
  conda: SuperCellMultiomics_pyenv
  benchmark:
    "benchmark/multiomicVelocity/{velocitySmp}/{runMv}/multivelo_benchmark.txt"
  shell: "python3 snakemakeWorkflows/multiomicVelocity/python/MultiveloCL.py -r {input.rna} -a {input.atac}\
          -w output/multiomicVelocity/{wildcards.velocitySmp}/{wildcards.runMv}/\
          -o output/multiomicVelocity/{wildcards.velocitySmp}/{wildcards.runMv}"


rule singleCell_multivelo:
  input: rna =  "input/multiomicVelocity/{velocitySmp}/adata_postpro_rna.h5ad",
         atac = "input/multiomicVelocity/{velocitySmp}/adata_postpro_atac.h5ad",
  output: "output/multiomicVelocity/{velocitySmp}/singleCell/multivelo_result.h5ad"
  conda: SuperCellMultiomics_pyenv
  params: adjLikelihood =  lambda w: config[w.velocitySmp]['adjLikelihood']
  benchmark:
    "benchmark/multiomicVelocity/{velocitySmp}/singleCell/multivelo_benchmark.txt"
  shell: "python3 snakemakeWorkflows/multiomicVelocity/python/MultiveloCL.py -r {input.rna} -a {input.atac}\
          -o output/multiomicVelocity/{wildcards.velocitySmp}/singleCell {params.adjLikelihood}"

