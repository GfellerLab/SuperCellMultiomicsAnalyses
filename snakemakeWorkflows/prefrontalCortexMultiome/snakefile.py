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
        outdir = OUTPUT_DIR + "/{n_smps}_samples/metacells/",
        python = config["prefrontalCortexMultiome"]["python"],
    singularity:
        config["prefrontalCortexMultiome"]["sif"]
    output: 
        OUTPUT_DIR + "/{n_smps}_samples/metacells/multiome_membership.qs2",
        OUTPUT_DIR + "/{n_smps}_samples/metacells/multiome_sample_ids.qs2",
        OUTPUT_DIR + "/{n_smps}_samples/metacells/multiome_keep_peaks.qs2"
    benchmark:
        OUTPUT_DIR + "/{n_smps}_samples/metacells/run_multiome_supercell2.txt"
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

rule run_multiome_integration_RNA_STACAS_unsup:
    input:
        sample_ids = OUTPUT_DIR + "/{n_smps}_samples/metacells/multiome_sample_ids.qs2",
        script = "snakemakeWorkflows/prefrontalCortexMultiome/R/run_multiome_integration_RNA_STACAS.R",
        functions = "snakemakeWorkflows/prefrontalCortexMultiome/R/functions/tools.r",
        integration_functions = "snakemakeWorkflows/prefrontalCortexMultiome/R/integration_method_functions.r"
    params:
        outdir = OUTPUT_DIR + "/{n_smps}_samples/RNA/",
        metacell_dir = OUTPUT_DIR + "/metacells_{n_smps}_samples/sample_cache" 
    singularity:
        config["prefrontalCortexMultiome"]["sif"]
    output:
        OUTPUT_DIR + "/{n_smps}_samples/RNA/multiome_RNA_integrated.qs2"
    benchmark:
        OUTPUT_DIR + "/{n_smps}_samples/RNA/run_multiome_integration_RNA_STACAS.txt"
    shell:
        """
        Rscript {input.script}  \
           --input-metacell-dir {params.metacell_dir} \
           --output-dir {params.outdir} --samples {input.sample_ids} \
           --n-core-samples 10 \
           --cache-threads 6 \
           --functions {input.functions} \
           --integration-r-fun {input.integration_functions}
       """


rule run_multiome_integration_ATAC_SignacRLSI_unsup:
    input:
        sample_ids = OUTPUT_DIR + "/{n_smps}_samples/metacells/multiome_sample_ids.qs2",
        script = "snakemakeWorkflows/prefrontalCortexMultiome/R/run_multiome_integration_ATAC_SignacRLSI.R",
        functions = "snakemakeWorkflows/prefrontalCortexMultiome/R/functions/tools.r",
        peaks = OUTPUT_DIR + "/{n_smps}_samples/metacells/multiome_keep_peaks.qs2"
    params:
        outdir = OUTPUT_DIR + "/{n_smps}_samples/ATAC/",
        metacell_dir = OUTPUT_DIR + "/{n_smps}_samples/metacells/sample_cache"
    singularity:
        config["prefrontalCortexMultiome"]["sif"]
    output:
        OUTPUT_DIR + "/{n_smps}_samples/ATAC/multiome_ATAC_integrated.qs2"
    benchmark:
        OUTPUT_DIR + "/{n_smps}_samples/ATAC/run_multiome_integration_ATAC_SignacRLSI.txt"
    shell:
        """
        Rscript {input.script}  \
           --input-metacell-dir {params.metacell_dir} \
           --keep_peaks {input.peaks} \
           --output-dir {params.outdir} --samples {input.sample_ids} \
           --n-core-samples 10 \
           --cache-threads 6 \
           --functions {input.functions} 
       """


rule run_multiome_integration_WNN_final_unsup:
    input:
        rna_integrated = OUTPUT_DIR + "/{n_smps}_samples/RNA/multiome_RNA_integrated.qs2",
        atac_integrated = OUTPUT_DIR + "/{n_smps}_samples/ATAC/multiome_ATAC_integrated.qs2",
        script = "snakemakeWorkflows/prefrontalCortexMultiome/R/run_multiome_integration_WNN_final.R",
        functions = "snakemakeWorkflows/prefrontalCortexMultiome/R/functions/tools.r"
    params:
        outdir = OUTPUT_DIR + "/{n_smps}_samples/WNN/",
        rna_reduc = "unsup_STACAS",
        atac_reduc= "unsup_SignacRLSI"
    singularity:
        config["prefrontalCortexMultiome"]["sif"]
    output:
        OUTPUT_DIR + "/{n_smps}_samples/WNN/multiome_integrated.qs2"
    benchmark:
        OUTPUT_DIR + "/{n_smps}_samples/WNN/run_multiome_integration_WNN_final.txt"
    shell:
       """
       Rscript {input.script} \
          --rna-integrated {input.rna_integrated} \
          --atac-integrated {input.atac_integrated} \
          --rna-reduc {params.rna_reduc} \
          --atac-reduc {params.atac_reduc} \
          --output-dir {params.outdir} 
       """

rule report:
    input:  expand(OUTPUT_DIR+"/{n_smps}_samples/WNN/multiome_integrated.qs2",n_smps = [10])
    output: "reports/prefrontalCortexMultiome/prefrontalCortexMultiome_analysis.html"
    shell:
        """
        touch {output}
        """





