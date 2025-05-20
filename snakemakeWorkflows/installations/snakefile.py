#############################################
######   package installations      #########
#############################################

# rule install_seuratData_package:
#   output: "config/installed_seuratData"
#   conda: SuperCellMultiomics
#   shell: "Rscript -e 'devtools::install_github(\"satijalab/seurat-data\")';\
#           touch {output}"
# 
# rule install_scDblFinder:
#     output: "config/scDblFinder_installed"
#     conda: SuperCellMultiomics
#     shell: "Rscript -e 'BiocManager::install(\"scDblFinder\")' && touch {output}"
#   
# rule install_ATAC_analysis_packages:
#   input: "config/SuperCellMultiomics/installedSuperCellMultiomics"
#   output: "config/installedAtacPackages"
#   conda: SuperCellMultiomics
#   shell : "Rscript -e 'remotes::install_github(\"satijalab/seurat-wrappers\");\
#                        remotes::install_github(\"mojaveazure/seurat-disk\");\
#                        BiocManager::install(\"chromVAR\");\
#                        BiocManager::install(\"JASPAR2020\");\
#                        BiocManager::install(\"TFBSTools\");\
#                        BiocManager::install(\"TFBSTools\");\
#                        BiocManager::install(\"motifmatchr\");\
#                        BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\");\
#                        BiocManager::install(\"biovizBase\");\
#                        setRepositories(ind=1:3);\
#                        install.packages(\"Signac\",repos = \"https://cloud.r-project.org\");\
#                        BiocManager::install(\"EnsDb.Hsapiens.v86\")';\
#                        touch  {output}"
#                        
# # remotes::install_github(\"rstudio/reticulate\",upgrade = \"never\"); #temporary fix for reading sparse matrix with R anndata https://github.com/rstudio/reticulate/issues/1417 \
# 
# 
# rule installRPackages:
#   output: "config/SuperCellMultiomics/installedPackages"
#   conda: SuperCellMultiomics
#   shell: "Rscript -e 'install.packages(c(\"igraph\", \"RANN\", \"WeightedCluster\", \"corpcor\", \"weights\", \"Hmisc\", \"Matrix\", \"matrixStats\", \"plyr\", \"irlba\", \"patchwork\", \"gtools\", \"ggplot2\", \"umap\", \"entropy\", \"Rtsne\", \"dbscan\", \"cowplot\", \"scales\", \"plotfunctions\", \"proxy\", \"Seurat\", \"data.table\", \"foreach\", \"doParallel\"),repos = \"https://cloud.r-project.org\")';\
#           cd config/SuperCellMultiomics;\
#           Rscript -e 'if (!require(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\",repos = \"https://cloud.r-project.org\");install.packages(\"ggsignif\",repos = \"https://cloud.r-project.org\"); \
#           remotes::install_github(\"carmonalab/STACAS\");\
#           BiocManager::install(\"bluster\")';\
#           cd ../../; touch {output}"
# 
# rule installSuperCellMultiomics:
#   input: "config/SuperCellMultiomics/installedPackages"
#   output: "config/SuperCellMultiomics/installedSuperCellMultiomics"
#   conda: SuperCellMultiomics
#   shell: "cd config/SuperCellMultiomics/;\
#           Rscript -e 'devtools::build();devtools::install_local(\"../SuperCellMultiomics_1.0.tar.gz\")';\
#           cd ../../;touch {output}"
# 
# rule install_MATK:
#   input: "config/SuperCellMultiomics/installedPackages"
#   output: "config/SuperCellMultiomics/installedMATK"
#   conda: SuperCellMultiomics
#   shell: "cd config/SuperCellMultiomics/;\
#           Rscript -e 'remotes::install_github(\"GfellerLab/SuperCell\",upgrade = \"never\");\
#           remotes::install_github(\"GfellerLab/MetacellAnalysisToolkit\",upgrade = \"never\");\
#           cd ../../;touch {output}"
#           
# rule install_bgzip:
#   output: "config/installed_bgzip"
#   conda: SuperCellMultiomics
#   params: workdir = wdir
#   shell: "cd config;wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2; \
#           tar -xf htslib-1.16.tar.bz2; \
#           cd htslib-1.16; \
#           ./configure --prefix={params.workdir}/config/bin/htslib-1.16; \
#           make; \
#           make install;\
#           touch installed_bgzip"

#rule create_jaspar_2024_mouse_pfm:
#  output: "input/motifs/mouse/jaspar2024_pfm.rds"
#  conda: JASPAR2024
#  shell: "Rscript -e 'library(JASPAR2024);\
#                      JASPAR2024 <- JASPAR2024();\
#                      db(JASPAR2024);\
#                      pfm <- TFBSTools::getMatrixSet(x = db(JASPAR2024),\
#                      opts = list(collection = \"CORE\", species = \"Mus musculus\",  matrixtype = \"PFM\", all_versions = FALSE));\
#                      id.version.df <- stringr::str_split_fixed(names(pfm),pattern=\"\\\\.\",n=2); \
#                      unique.id <- c(); \
#                      for (i in unique(id.version.df[,1])) {{ \
#                        latest <- max(id.version.df[id.version.df[,1]==i,2]); \
#                        latest.id <- paste0(i,\".\",latest); \
#                        unique.id <- c(unique.id,latest.id); \
#                      }}; \
#                      pfm <- pfm[unique.id]; \
#                      saveRDS(pfm,\"{output}\")'"
                      
                      
#rule create_jaspar_2024_human_pfm:
#  output: "input/motifs/human/jaspar2024_pfm.rds"
#  conda: JASPAR2024
#  shell: "Rscript -e 'library(JASPAR2024);\
#                      JASPAR2024 <- JASPAR2024();\
#                      db(JASPAR2024);\
#                      pfm <- TFBSTools::getMatrixSet(x = db(JASPAR2024),\
#                      opts = list(collection = \"CORE\", species = 9606,  matrixtype = \"PFM\", all_versions = TRUE));\
#                      id.version.df <- stringr::str_split_fixed(names(pfm),pattern=\"\\\\.\",n=2); \
#                      unique.id <- c(); \
#                      for (i in unique(id.version.df[,1])) {{latest <- max(id.version.df[id.version.df[,1]==i,2]);latest.id <- paste0(i,\".\",latest); unique.id <- c(unique.id,latest.id)}}; \
#                      pfm <- pfm[unique.id]; \
#                      saveRDS(pfm,\"{output}\")'"
                      
# rule install_sup_tools:
#   output: "config/installed_sup_tools"
#   conda: SupervisedTools
#   shell: "Rscript -e 'if (!require(\"BiocManager\", quietly = TRUE)){{install.packages(\"BiocManager\",repos = \"https://cloud.r-project.org\")}};\
#                       BiocManager::install(\"UCell\",update = F);\
#                       install.packages(\"scGate\",repos = \"https://cloud.r-project.org\");\
#                       if(!require(dlm)) devtools::install_version(\"dlm\", version = \"1.1.5\", repos = \"http://cran.us.r-project.org\",upgrade = \"never\"); \
#                       if(!require(Rmagic)) devtools::install_version(\"Rmagic\", version = \"2.0.3\", repos = \"http://cran.us.r-project.org\",upgrade = \"never\"); \
#                       if(!require(cutoff.scATOMIC)) devtools::install_github(\"inofechm/cutoff.scATOMIC\", force = T,upgrade = \"never\");\
#                       if(!require(scATOMIC)) devtools::install_github(\"abelson-lab/scATOMIC\",upgrade = \"never\");\
#                       library(reticulate);\
#                       library(Rmagic);\
#                       install.magic();\
#                       pymagic_is_available()';\
#           Rscript -e 'install.packages(c(\"igraph\", \"RANN\", \"WeightedCluster\", \"corpcor\", \"weights\", \"Hmisc\", \"Matrix\", \"matrixStats\", \"plyr\", \"irlba\", \"patchwork\", \"gtools\", \"ggplot2\", \"umap\", \"entropy\", \"Rtsne\", \"dbscan\", \"cowplot\", \"scales\", \"plotfunctions\", \"proxy\", \"Seurat\", \"data.table\", \"foreach\", \"doParallel\"),repos = \"https://cloud.r-project.org\")';\
#           cd config/SuperCellMultiomics;\
#           Rscript -e 'install.packages(\"ggsignif\",repos = \"https://cloud.r-project.org\"); \
#           remotes::install_github(\"carmonalab/STACAS\",upgrade = \"never\");\
#           BiocManager::install(\"bluster\",update = F)';\
#           Rscript -e 'devtools::build();devtools::install_local(\"../SuperCellMultiomics_1.0.tar.gz\")';\
#           cd ../../;touch {output}"
                      
# rule create_jaspar_2024_vertebrates_pfm:
#   output: "input/JASPAR2024_motifs/vertebrate_pfm.rds"
#   conda: JASPAR2024
#   shell: "Rscript -e 'library(JASPAR2024);\
#                       JASPAR2024 <- JASPAR2024()\
#                       db(JASPAR2024);\
#                       pfm <- TFBSTools::getMatrixSet(x = db(JASPAR2024),\
#                       opts = list(collection = \"CORE\", tax_group = \"vertebrates\",  matrixtype = \"PFM\", all_versions = FALSE))\
#                       saveRDS(pfm,\"{output}\")'"

# rule install_seurat_dataset:
#   input: "config/installed_seuratData","config/installedAtacPackages"
#   output: "config/installed_{seuratDataset}"
#   conda: SuperCellMultiomics
#   shell: "Rscript -e 'SeuratData::InstallData(\"{wildcards.seuratDataset}\")'; touch {output}"
#   
# rule install_seurat_v5:
#   output: "config/installed_seuratV5"
#   conda: seuratV5
#   shell: "Rscript -e 'remotes::install_github(\"satijalab/seurat\", \"seurat5\", quiet = TRUE);\
#                       remotes::install_github(\"bnprks/BPCells\", quiet = TRUE);\
#                       remotes::install_github(\"carmonalab/STACAS\");\
#                       remotes::install_github(\"satijalab/seurat-wrappers\", \"seurat5\", quiet = TRUE)' && touch {output}"            
