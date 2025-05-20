#### Postate cancer analysis

prostateSamples = list(config["prostateHanSamples"])
cleanProstateSamples = list(config["cleanProstateHanSamples"])

pathCellRangerProstate = config["pathCellRangerProstateData"]

#workdir: "../../"
  
rule RNA_object_construction:
  input: "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_3/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_4/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/6M/filtered_feature_bc_matrix.h5"

  output:"input/prostateCancer10xMultiome/rna.All.combined.integrated.Rdata"
  benchmark: "benchmark/10xMultiomeProstateHan22/single_cell/RNA_object_construction.txt"
  conda: SuperCellMultiomics
  shell: "cd input/prostateCancer10xMultiome;Rscript ../../snakemakeWorkflows/10xMultiomeProstateHan22/R/RNA_object_construction_prostateCancer.R"

rule ATAC_object_construction:
  input: "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_3/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_4/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_1/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_2/filtered_feature_bc_matrix.h5",
         "input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/6M/filtered_feature_bc_matrix.h5"

  output:"input/prostateCancer10xMultiome/atac.prostate.merge.Rdata"
  conda: SuperCellMultiomics
  benchmark: "benchmark/10xMultiomeProstateHan22/single_cell/ATAC_object_construction.txt"
  shell: "cd input/prostateCancer10xMultiome;Rscript ../../snakemakeWorkflows/10xMultiomeProstateHan22/R/ATAC_object_construction_prostateCancer.R"

rule ATAC_RNA_integration:
  input: "input/prostateCancer10xMultiome/rna.All.combined.integrated.Rdata","input/prostateCancer10xMultiome/atac.prostate.merge.Rdata"
  output : "input/prostateCancer10xMultiome/wnn.all.combined.integrated_meta.data.csv"
  conda: SuperCellMultiomics
  benchmark: "benchmark/10xMultiomeProstateHan22/single_cell/ATAC_RNA_integration.txt"
  shell: "cd input/prostateCancer10xMultiome;Rscript ../../snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_RNA_ATAC_prostateCancer.R"


############################################
######## Analysis with SuperCell ###########
############################################

rule GetConsensusPeaks:
  conda: SuperCellMultiomics
  input: expand("input/prostateCancer10xMultiome/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/{sample}/filtered_feature_bc_matrix.h5",sample = prostateSamples)
  output: consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/GetConsensusPeaks_prostateCancer.R"
  
rule call_macs2_peaks_smp:
  input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
         index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
         h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
         packages=  "config/installedAtacPackages"
  output: "output/10xMultiomeProstateHan22/macs2_peaks/{sample}/macs2_peaks.rds"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/call_macs2_peaks_CL.R\
          -i {input.h5}\
          -f {input.fragment}\
          -o output/10xMultiomeProstateHan22/macs2_peaks/{wildcards.sample}"
  

rule GetConsensusPeaks_macs2:
  conda: SuperCellMultiomics
  input: expand("output/10xMultiomeProstateHan22/macs2_peaks/{sample}/macs2_peaks.rds",sample = cleanProstateSamples)
  output: consensusPeaks = "output/10xMultiomeProstateHan22/macs2_peaks/peakset_macs2.csv"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/GetConsensusPeaks_macs2_prostateCancer.R\
          -i output/10xMultiomeProstateHan22/macs2_peaks/\
          -o output/10xMultiomeProstateHan22/macs2_peaks/"
  
rule metacellPerSample:
  input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
         index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
         h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
         packages=  "config/installedAtacPackages",
         consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
  output: "output/10xMultiomeProstateHan22/logNorm/{sample}/metacells.rds"
  params: python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/metacellPerSampleCL.R -i {input.h5} \
         -f {input.fragment} \
         -o output/10xMultiomeProstateHan22/logNorm/{wildcards.sample}\
         -p {params.python}\
         -c {input.consensusPeaks} \
         -g 10"
         
# tested not retained  
rule metacellPerSample_macs2_peak_set:
  input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
         index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
         h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
         packages=  "config/installedAtacPackages",
         macsPeak = "output/10xMultiomeProstateHan22/logNorm/macs2_major_type/macs2.peaks.rds"
  output: "output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/{sample}/metacells.rds"
  params: python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/metacellPerSampleCL.R -i {input.h5} \
         -f {input.fragment} \
         -o output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/{wildcards.sample}\
         -p {params.python}\
         -a {input.macsPeak} \
         -g 10"
         
# rule metacellPerSample_macs2_peak:
#   input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
#          index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
#          h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
#          packages=  "config/installedAtacPackages",
#          consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
#   output: "output/10xMultiomeProstateHan22/logNorm/macs2/{sample}/metacells.rds"
#   params: python = Seacells + "/bin/python"
#   conda: SuperCellMultiomics
#   shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/metacellPerSampleCL.R -i {input.h5} \
#          -f {input.fragment} \
#          -o output/10xMultiomeProstateHan22/logNorm/macs2/{wildcards.sample}\
#          -p {params.python}\
#          -c {input.consensusPeaks} \
#          -g 10 \
#          -m .snakemake/conda/29a8d0b0/bin/macs2"
         
# tested not retained  
rule metacellPerSample_g5:
  input: fragment = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz",
         index = pathCellRangerProstate + "{sample}/atac_fragments.tsv.gz.tbi",
         h5 = pathCellRangerProstate + "{sample}/filtered_feature_bc_matrix.h5",
         packages=  "config/installedAtacPackages",
         consensusPeaks = "input/prostateCancer10xMultiome/peakset_minOverlap_0.5_remove.csv"
  output: "output/10xMultiomeProstateHan22/gamma_5/logNorm/{sample}/metacells.rds"
  params: python = Seacells + "/bin/python"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/metacellPerSampleCL.R -i {input.h5} \
         -f {input.fragment} \
         -o output/10xMultiomeProstateHan22/gamma_5/logNorm/{wildcards.sample}\
         -p {params.python}\
         -c {input.consensusPeaks} \
         -g 5"

# tested not retained  
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
  input: expand("output/10xMultiomeProstateHan22/logNorm/{sample}/metacells.rds",sample = prostateSamples)
  output: "output/10xMultiomeProstateHan22/logNorm/combined.metacells.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/10xMultiomeProstateHan22/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R -o output/10xMultiomeProstateHan22/logNorm/ -i harmony"


rule annotate_metacell_cluster:
  input: "output/10xMultiomeProstateHan22/logNorm/combined.metacells.rds"
  output: immuneMetacells = "output/10xMultiomeProstateHan22/logNorm/immune.combined.metacells.rds"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Annotating_metacells.R -i {input} -o output/10xMultiomeProstateHan22/logNorm/"

rule peak_calling_immune_chromvar_jaspar_2024_mouse:
  input: seurat = "output/10xMultiomeProstateHan22/logNorm/immune.combined.metacells.rds",
         anno = "output/10xMultiomeProstateHan22/logNorm/combined.metacells_all_anno.csv",
         motifs = "input/motifs/mouse/jaspar2024_pfm.rds"
  output: "output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune_JASPAR_motifs/combined.metacells.with.macs.peak.rds"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/macs_peak_calling.R -i {input.seurat}\
  -m {input.anno} \
  -o output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune_JASPAR_motifs/\
  -p {input.motifs}"

rule peak_calling_immune_chromvar_hoccomoco:
  input: seurat = "output/10xMultiomeProstateHan22/logNorm/immune.combined.metacells.rds",
         anno = "output/10xMultiomeProstateHan22/logNorm/combined.metacells_all_anno.csv",
         motifs = "input/HOCCOMOCO_motifs/H12CORE_pfm.rds"
  output: "output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune_hoccomoco_motifs/combined.metacells.with.macs.peak.rds"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/macs_peak_calling.R -i {input.seurat}\
  -m {input.anno} \
  -o output/10xMultiomeProstateHan22/logNorm/macs2_major_type/immune_hoccomoco_motifs/\
  -p {input.motifs}"


rule peak_calling_all_major_type:
  input: seurat = "output/10xMultiomeProstateHan22/logNorm/combined.metacells.rds",
         anno = "output/10xMultiomeProstateHan22/logNorm/combined.metacells_all_anno.csv"
  output: "output/10xMultiomeProstateHan22/logNorm/macs2_major_type/macs2.peaks.rds"
  conda: SuperCellMultiomics
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/macs_peak_calling.R -i {input.seurat} -m {input.anno} -u peaks -o output/10xMultiomeProstateHan22/logNorm/macs2_major_type/"


  
rule Integrating_metacells_RNA_ATAC_prostateCancer_macs2_peak_set:
  input: expand("output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/{sample}/metacells.rds",sample = cleanProstateSamples)
  output: "output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/combined.metacells.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/10xMultiomeProstateHan22/macs2_peak_set/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R\
  -o output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/\
  -i harmony"
  
rule annotate_metacell_cluster_metacell_macs2_peak_set: 
  input: sobj = "output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/combined.metacells.rds",
         rmd = "reports/10xMultiomeProstateHan22/results_metacell_macs_peaks_major_type.Rmd"
  output: "output/10xMultiomeProstateHan22/logNorm/macs2_peak_set/immune.combined.metacells.rds"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'rmarkdown::render(\"{input.rmd}\")'"


# rule Integrating_metacells_RNA_ATAC_prostateCancer_macs2:
#   input: expand("output/10xMultiomeProstateHan22/logNorm/macs2/{sample}/metacells.rds",sample = cleanProstateSamples)
#   output: "output/10xMultiomeProstateHan22/logNorm/macs2/combined.metacells.rds"
#   conda: SuperCellMultiomics
#   benchmark:
#     "benchmark/10xMultiomeProstateHan22/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
#   shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R -o output/10xMultiomeProstateHan22/logNorm/macs2/ -i harmony"

# tested not retained  
rule Integrating_metacells_RNA_ATAC_prostateCancer_gamma_5:
  input: expand("output/10xMultiomeProstateHan22/gamma_5/logNorm/{sample}/metacells.rds",sample = prostateSamples)
  output: "output/10xMultiomeProstateHan22/gamma_5/logNorm/combined.metacells.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/10xMultiomeProstateHan22/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R -o output/10xMultiomeProstateHan22/gamma_5/logNorm/ -i harmony"


# tested not retained  
rule Integrating_metacells_RNA_ATAC_prostateCancer_SCT:
  input: expand("output/10xMultiomeProstateHan22/SCT/{sample}/metacells.rds",sample = prostateSamples)
  output: "output/10xMultiomeProstateHan22/SCT/combined.metacells.rds"
  conda: SuperCellMultiomics
  benchmark:
    "benchmark/10xMultiomeProstateHan22/SCT/Integrating_metacells_RNA_ATAC_prostateCancer.txt"
  shell: "Rscript snakemakeWorkflows/10xMultiomeProstateHan22/R/Integrating_metacells_RNA_ATAC_prostateCancer.R -o output/10xMultiomeProstateHan22/SCT/ -n SCT -f 0.4"
