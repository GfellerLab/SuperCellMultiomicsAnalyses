#### Postate cancer analysis

prostateSamples = list(config["prostateHanSamples"])
pathCellRangerProstate = config["pathCellRangerProstateData"]

#workdir: "../../"
  
# rule RNA_object_construction:
#   input: "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_3/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_4/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/6M/filtered_feature_bc_matrix.h5"
#          
#   output:"input/prostateCancer10xMultiome/rna.All.combined.integrated.Rdata"
#   conda: han_env
#   shell: "Rscript R/RNA_object_construction_prostateCancer.R"
#   
# rule ATAC_object_construction:
#   input: "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_3/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_4/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_1/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_2/filtered_feature_bc_matrix.h5",
#          "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/6M/filtered_feature_bc_matrix.h5"
#          
#   output:"input/prostateCancer10xMultiome/atac.prostate.merge.Rdata"
#   conda: SuperCellMultiomics
#   shell: "Rscript R/ATAC_object_construction_prostateCancer.R"
#   
# rule ATAC_RNA_integration:
#   input: "input/prostateCancer10xMultiome/rna.All.combined.integrated.Rdata","input/prostateCancer10xMultiome/atac.prostate.merge.Rdata"
#   output : "input/prostateCancer10xMultiome/wnn.all.combined.integrated_meta.data.csv"
#   conda: SuperCellMultiomics
#   shell: "Rscript R/Integrating_RNA_ATAC_prostateCancer.R"


############################################
######## Analysis with SuperCell ###########
############################################

rule GetConsensusPeaks:
  conda: SuperCellMultiomics
  input: expand("input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/{sample}/filtered_feature_bc_matrix.h5",sample = prostateSamples)
  output: consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
  shell: "Rscript R/GetConsensusPeaks_prostateCancer.R"
  
rule metacellPerSample:
  input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
         index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
         h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
         packages=  "config/installedAtacPackages",
         consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
  output: "output/10xMultiomeProstateHan22/{sample}/metacells.rds"
  params: python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/metacellPerSampleCL.R -i {input.h5} \
         -f {input.fragment} \
         -o output/10xMultiomeProstateHan22/{wildcards.sample}\
         -p {params.python}\
         -c {input.consensusPeaks} \
         -g 10"

rule metacellPerSample_SCT:
  input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
         index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
         h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
         packages=  "config/installedAtacPackages",
         consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
  output: "output/10xMultiomeProstateHan22/SCT/{sample}/metacells.rds"
  params: python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/metacellPerSampleCL.R -i {input.h5} \
         -f {input.fragment} \
         -o output/10xMultiomeProstateHan22/SCT/{wildcards.sample}\
         -p {params.python}\
         -n SCT \
         -c {input.consensusPeaks} \
         -g 10"
         
rule Integrating_metacells_RNA_ATAC_prostateCancer:
  input: expand("output/10xMultiomeProstateHan22/{sample}/metacells.rds",sample = prostateSamples)
  output: "output/10xMultiomeProstateHan22/combined.metacells.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/10xMultiomeProstateHan22/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R"
  
rule Integrating_metacells_RNA_ATAC_prostateCancer_SCT:
  input: expand("output/10xMultiomeProstateHan22/SCT/{sample}/metacells.rds",sample = prostateSamples)
  output: "output/10xMultiomeProstateHan22/SCT/combined.metacells.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/10xMultiomeProstateHan22/SCT/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R -o output/10xMultiomeProstateHan22/SCT/ -n SCT -f 0.4"
