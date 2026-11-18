
rule download_data:
    output: 
        rna_h5ad = INPUT_DIR + "/final_rna_data.h5ad",
        atac_h5ad =  INPUT_DIR + "/final_rna_data.h5ad"
    params: 
        rna_url = config["prefrontalCortexMultiome"]["rna_url"],
        atac_url = config["prefrontalCortexMultiome"]["atac_url"]
    shell: 
        "cd {params.input_dir} && wget {params.rna_url} && wget {params.atac_url}"
        
rule run_multiome_supercell2_unssup:
    input: 
        script = "snakemakeWorkflows/prefrontalCortexMultiome/R/run_multiome_supercell2.R"
        functions = "snakemakeWorkflows/prefrontalCortexMultiome/R/functions/run_multiome_supercell2.R"
        rna_h5ad = INPUT_DIR + "/final_rna_data.h5ad",
        atac_h5ad =  INPUT_DIR + "/final_atac_data.h5ad"
    params: 
        outdir = OUTPUT_DIR + "/metacells_{n_smps}_samples/"
        python = config["prefrontalCortexMultiome"]["python"]
        r_fun = config["prefrontalCortexMultiome"]["r_fun"]
    singularity:
        config["prefrontalCortexMultiome"]["supercell_qs2_sif"]
    output: 
        OUTPUT_DIR + "/metacells_{n_smps}_samples/multiome_membership.qs2"
    benchmark:
        OUTPUT_DIR + "/metacells_{n_smps}_samples/run_multiome_supercell2.txt"
    shell:
        """
        Rscript  \
            --rna-h5ad {input.rna_h5ad} --atac-h5ad {input.atac_h5ad} \
            --output-dir {params.outdir} \
            --python {params.python} \
            --functions {input.functions}
            --n-samples {wildcards.samples} \
            --min-peak-fraction 0.05 \
            --cache-threads 1 
        """

rule report:
    input:  expand(OUTDIR+"/metacells_{n_smps}_samples/multiome_membership.qs2",n_smps = [10,50])
    output: "reports/prefrontalCortexMultiome/prefrontalCortexMultiome_analysis.html"
    shell:
        """
        touch {output}
        """





