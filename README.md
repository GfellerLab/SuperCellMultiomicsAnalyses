# Analysis Workflows for our study : SuperCell2.0 enables semi-supervised construction of multimodal metacell atlases

* Benchmark Supercell2.0 against [SEACells](https://github.com/dpeerlab/SEACells) and [metacell-2](https://www.weizmann.ac.il/math/tanay/research-activities/metacell-2) on single-cell multimodal data. 
* Analyse intermodality consistency at the metacell level.
* Benchmark the construction of the PBMC CITE-seq metacell atlas with semi-supervised and unsupervised workflows using our metacell tool [SuperCell2.0](https://github.com/GfellerLab/SuperCell) and [STACAS](https://github.com/carmonalab/STACAS).
* Annotate and analyze the semi-supervised PBMC CITE-seq metacell atlas.

## Additional data

Additional data (eg. containers, transcription factor motifs pfm matrix, gene signatures from the literature) required for this workflow are availbale on [zenodo](https://doi.org/10.5281/zenodo.18613560)

## Launch the snakemake workflow (snakemake v7.15.2) on a cluster with slurm and singularity

### Workflow env

This workflow has been developped using snakemake v7.15.2.
You can install our workflow conda env like this:

    conda env create -n snakemake -f config/snakemake_v7.15.2.yml

#### load singularity module

    module load singularityce
   

#### launch snakemake revisions
    
    conda activate snakemake
    
    snakemake -j 30 \
        -kps snakefile_manuscript.py \
        --configfile config/workflow_revisions.yml \
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


#### launch snakemake
    
    conda activate snakemake
    
    snakemake -j 30 \
        -kps snakefile_manuscript.py \
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

