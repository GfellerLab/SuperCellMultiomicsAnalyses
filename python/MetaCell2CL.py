import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import metacells as mc
import numpy as np
import os
import pandas as pd
import scipy.sparse as sp
import statistics
from math import hypot
from scipy import sparse
import sys, getopt
from pathlib import Path
np.random.seed(123456)

def main(argv):
    matrixFile = ''
    outDir = ''
    try:
        opts, args = getopt.getopt(argv,"i:g:d:r:o:",["inputH5ad=","gamma=","outDir="])
    except getopt.GetoptError:
        print('MetaCell2CL.py -i <inputH5ad> -g <gamma> -o <outDir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('MetaCell2CL.py -i <inputH5ad> -g <gamma> -o <outDir>')
            sys.exit()
        elif opt in ("-i", "--inputH5ad"):
            inputH5ad = arg
        elif opt in ("-g", "--gamma"):
            gamma = float(arg)
        elif opt in ("-o", "--outDir"):
            print(arg)
            outDir = arg

    print('Output dir is "', outDir)
    print('inputH5ad is "', inputH5ad)
    print('gamma is "', gamma)
    adata = sc.read(inputH5ad)
    raw  = ad.AnnData(X = adata.raw.X,
                  obs = adata.obs,
                  var = adata.raw.var)
    proj_name = "MetaCell2RNA"
    # set MC object/project name
    mc.ut.set_name(raw, proj_name)
    print(raw.shape)
    excluded_gene_names = [] # for example, ['IGHMBP2', 'IGLL1', 'IGLL5', 'IGLON5', 'NEAT1', 'TMSB10', 'TMSB4X']
    excluded_gene_patterns = ['MT-.*']
    ### The first round (high/low UMIs), very lenient thresholds, we assume the dataset has already been filtered (annotation for each cell)
    properly_sampled_min_cell_total = 0
    properly_sampled_max_cell_total = pow(10,9)
    total_umis_of_cells = mc.ut.get_o_numpy(raw, name='__x__', sum=True)

    too_small_cells_count = sum(total_umis_of_cells < properly_sampled_min_cell_total)
    too_large_cells_count = sum(total_umis_of_cells > properly_sampled_max_cell_total)

    too_small_cells_percent = 100.0 * too_small_cells_count / len(total_umis_of_cells)
    too_large_cells_percent = 100.0 * too_large_cells_count / len(total_umis_of_cells)

    print(f"Will exclude %s (%.2f%%) cells with less than %s UMIs"
      % (too_small_cells_count,
         too_small_cells_percent,
         properly_sampled_min_cell_total))
    print(f"Will exclude %s (%.2f%%) cells with more than %s UMIs"
      % (too_large_cells_count,
         too_large_cells_percent,
         properly_sampled_max_cell_total))
    mc.pl.analyze_clean_genes(raw,
                          excluded_gene_names=excluded_gene_names,
                          excluded_gene_patterns=excluded_gene_patterns,
                          random_seed=123456)
    mc.pl.pick_clean_genes(raw)
    ## The second round (content of non-clean genes, e.g., mito-genes)
    properly_sampled_max_excluded_genes_fraction = 1 #we assume the dataset has already been filtered (annotation for each cell)


    excluded_genes_data = mc.tl.filter_data(raw, var_masks=['~clean_gene'])[0]
    excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
    excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells

    too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)

    too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)

    print(f"Will exclude %s (%.2f%%) cells with more than %.2f%% excluded gene UMIs"
      % (too_excluded_cells_count,
         too_excluded_cells_percent,
         100.0 * properly_sampled_max_excluded_genes_fraction))

    mc.pl.analyze_clean_cells(
        raw,
        properly_sampled_min_cell_total=properly_sampled_min_cell_total,
        properly_sampled_max_cell_total=properly_sampled_max_cell_total,
        properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction
    )
    mc.pl.pick_clean_cells(raw)
    # Extract clean dataset (with fillered cells and genes)
    clean = mc.pl.extract_clean_data(raw)

    # Estimate target_metacell_size(gamma):
    print(f'The requested graining level is {gamma}, lets estimate the target_metacell_size that should result in such graining level.')

    scale = 1 # incres or decrease if the obtained graining level (`gamma_obtained`) is significantly > or < then the requested one `gamma`

    N_c = clean.shape[0]

    # estimated mean UMI content in dowsampled data
    est_downsample_UMI = np.quantile(np.array(total_umis_of_cells), 0.5)

    target_metacell_size = int(est_downsample_UMI*gamma) * scale

    #target_metacell_size = int(gamma*np.mean(np.array(total_umis_of_cells))* scale)
    target_metacell_size

    mc.pl.divide_and_conquer_pipeline(
        clean,
        #feature_gene_names   = feature_gene_names, # comment this line to allow Metacell2 selecting features
        #forbidden_gene_names = forbidden_gene_names, # comment this line to allow Metacell2 selecting features
        target_metacell_size = target_metacell_size,
        random_seed = 123456)

    ## make anndata of metacells
    metacells = mc.pl.collect_metacells(clean, name='cell_lines.metacells')

    gamma_obtained = clean.shape[0]/metacells.shape[0]

    gamma_obtained

    clean.obs['membership'] = [np.nan for i in clean.obs.metacell]
    outlier = 0
    for m in range(0,len(clean.obs.metacell)):
        if clean.obs.metacell[m] >= 0:
            clean.obs['membership'][m] = clean.obs.metacell[m]+1
        else:
            clean.obs['membership'][m] = outlier-1 #outlier metacell of size one will receive negative membership
            outlier -=1

    clean.obs['membership'].to_csv(outDir+"/MetaCellMemberships.csv")

if __name__ == "__main__":
    main(sys.argv[1:])
