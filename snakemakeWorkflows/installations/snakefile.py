#############################################
######   package installations      #########
#############################################

rule install_seuratData_package:
  output: "config/installed_seuratData"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'devtools::install_github(\"satijalab/seurat-data\")';\
          touch {output}"
  
  
rule install_ATAC_analysis_packages:
  input: "config/SuperCellMultiomics/installedSuperCellMultiomics"
  output: "config/installedAtacPackages"
  conda: SuperCellMultiomics
  shell : "Rscript -e 'remotes::install_github(\"satijalab/seurat-wrappers\");\
                       remotes::install_github(\"mojaveazure/seurat-disk\");\
                       BiocManager::install(\"chromVAR\");\
                       BiocManager::install(\"JASPAR2020\");\
                       BiocManager::install(\"TFBSTools\");\
                       BiocManager::install(\"TFBSTools\");\
                       BiocManager::install(\"motifmatchr\");\
                       BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\");\
                       BiocManager::install(\"biovizBase\");\
                       setRepositories(ind=1:3);\
                       install.packages(\"Signac\",repos = \"https://cloud.r-project.org\");\
                       BiocManager::install(\"EnsDb.Hsapiens.v86\")';\
                       touch  {output}"
                       


rule installRPackages:
  output: "config/SuperCellMultiomics/installedPackages"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'install.packages(c(\"igraph\", \"RANN\", \"WeightedCluster\", \"corpcor\", \"weights\", \"Hmisc\", \"Matrix\", \"matrixStats\", \"plyr\", \"irlba\", \"patchwork\", \"gtools\", \"ggplot2\", \"umap\", \"entropy\", \"Rtsne\", \"dbscan\", \"cowplot\", \"scales\", \"plotfunctions\", \"proxy\", \"Seurat\", \"data.table\", \"foreach\", \"doParallel\"),repos = \"https://cloud.r-project.org\")';\
          cd config/SuperCellMultiomics;\
          Rscript -e 'if (!require(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\",repos = \"https://cloud.r-project.org\"); \
          BiocManager::install(\"bluster\")';\
          cd ../../; touch {output}"

rule installSuperCellMultiomics:
  input: "config/SuperCellMultiomics/installedPackages"
  output: "config/SuperCellMultiomics/installedSuperCellMultiomics"
  conda: SuperCellMultiomics
  shell: "cd config/SuperCellMultiomics/;\
          Rscript -e 'devtools::build();devtools::install_local(\"../SuperCellMultiomics_1.0.tar.gz\")';\
          cd ../../;touch {output}"
          
rule install_bgzip:
  output: "config/installed_bgzip"
  conda: SuperCellMultiomics
  params: workdir = wdir
  shell: "cd config;wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2; \
          tar -xf htslib-1.16.tar.bz2; \
          cd htslib-1.16; \
          ./configure --prefix={params.workdir}/config/bin/htslib-1.16; \
          make; \
          make install;\
          touch installed_bgzip"

                       
rule install_seurat_dataset:
  input: "config/installed_seuratData","config/installedAtacPackages"
  output: "config/installed_{seuratDataset}"
  conda: SuperCellMultiomics
  shell: "Rscript -e 'SeuratData::InstallData(\"{wildcards.seuratDataset}\")'; touch {output}"
