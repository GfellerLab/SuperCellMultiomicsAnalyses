# SuperCellMultiomicsAnalyses

Workflow for SuperCellMultiomics Analyses

## installation
conda install mamba -n base -c conda-forge

mkdir -p cluster/snakemake

cd config

git clone 


### Command line to test the workflow (dry run) 
    conda activate snakemake
    
    snakemake -j 1 \
        -nps snakefile.py \
        --configfile config/workflow.yml \
        --use-conda \ 
    
        
        
## Draw workflow 
    
    snakemake -j 1  -nps snakefile.py --configfile config/workflow.yml --use-conda --conda-frontend conda  --forceall --rulegraph > dag.dot
        
### Print workflow

    dot -Tpng dag.dot > dag.png

### command line to launch the snakemake workflow (snakemake v7.15.2) on a cluster with slurm and singularity

#### load singularity module

    module load singularityce
   

#### launch snakemake

   snakemake -j 30 \
        -kps snakefile.py \
        --configfile config/workflow.yml \
        --use-singularity \
        --conda-frontend conda \
        --cluster-config config/cluster.yml \
        --cluster "sbatch -A {cluster.account} \
        -p {cluster.partition} \
        -N {cluster.N} \
        -t {cluster.time} \
        --job-name {cluster.name} \
        --mem {cluster.mem} \
        --cpus-per-task {cluster.cpus-per-task}\
        --output {cluster.output} \
        --error {cluster.error}"

### command line to launch the snakemake workflow on a cluster with slurm

    conda activate snakemake

    export PATH="/dcsrsoft/spack/hetre/v1.2/spack/opt/spack/linux-rhel8-zen2/gcc-9.3.0/miniconda3-4.9.2-jle3zxexkucvzivrdnjax2kreaql6izj/bin:$PATH"
    snakemake -j 30 \
        -kps snakefile.py \
        --configfile config/workflow.yml \
        --use-conda \
        --conda-frontend conda \
        --cluster-config config/cluster.yml \
        --cluster "sbatch -A {cluster.account} \
        -p {cluster.partition} \
        -N {cluster.N} \
        -t {cluster.time} \
        --job-name {cluster.name} \
        --mem {cluster.mem} \
        --cpus-per-task {cluster.cpus-per-task}\
        --output {cluster.output} \
        --error {cluster.error}"
        

