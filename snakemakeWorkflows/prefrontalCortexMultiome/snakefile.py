INPUT_DIR = config["prefrontalCortexMultiome"]["input_dir"]
OUTPUT_DIR = config["prefrontalCortexMultiome"]["output_dir"]

rule download_data:
    output: 
        rna_h5ad = INPUT_DIR + "/final_rna_data.h5ad",
        atac_h5ad =  INPUT_DIR + "/final_atac_data.h5ad"
    params: 
        rna_url = config["prefrontalCortexMultiome"]["url_rna"],
        atac_url = config["prefrontalCortexMultiome"]["url_atac"]
    shell: 
        "cd {params.input_dir} && wget {params.rna_url} && wget {params.atac_url}"
        
rule run_multiome_supercell2_unsup:
    input: 
        script = "snakemakeWorkflows/prefrontalCortexMultiome/R/run_multiome_supercell2.R",
        functions = "snakemakeWorkflows/prefrontalCortexMultiome/R/functions/tools.r",
        rna_h5ad = INPUT_DIR + "/final_rna_data.h5ad",
        atac_h5ad =  INPUT_DIR + "/final_atac_data.h5ad"
    params: 
        outdir = OUTPUT_DIR + "/metacells_{n_smps}_samples/",
        python = config["prefrontalCortexMultiome"]["python"],
    singularity:
        config["prefrontalCortexMultiome"]["sif"]
    output: 
        OUTPUT_DIR + "/metacells_{n_smps}_samples/multiome_membership.qs2"
    benchmark:
        OUTPUT_DIR + "/metacells_{n_smps}_samples/run_multiome_supercell2.txt"
    shell:
        """
        Rscript {input.script}  \
            --rna-h5ad {input.rna_h5ad} --atac-h5ad {input.atac_h5ad} \
            --output-dir {params.outdir} \
            --python {params.python} \
            --functions {input.functions} \
            --n-samples {wildcards.n_smps} \
            --min-peak-fraction 0.05 \
            --cache-threads 1 
        """

rule report:
    input:  expand(OUTPUT_DIR+"/metacells_{n_smps}_samples/multiome_membership.qs2",n_smps = [10])
    output: "reports/prefrontalCortexMultiome/prefrontalCortexMultiome_analysis.html"
    shell:
        """
        touch {output}
        """





