import numpy as np
import pandas as pd
import scanpy as sc
import multivelo as mv
import scvelo as scv
import scanpy as sc
import os
import sys, getopt
from pathlib import Path

def main(argv):
    matrixFile = ''
    outDir = ''
    adjLikelihood = False
    try:
        opts, args = getopt.getopt(argv,"r:a:w:o:l",["rnaInputH5ad=","atacInputH5ad=","wnnDir=","outDir=","adjLikelihood"])
    except getopt.GetoptError:
        print('MultiveloCL.py -r <rnaInputH5ad> -a <atacInputH5ad> -w <wnnDir> -o <outDir> -l <adjLikelihood>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('MultiveloCL.py -r <rnaInputH5ad> -a <atacInputH5ad> -w <wnnDir> -l <adjLikelihood> -o <outDir> -l <adjLikelihood>')
            sys.exit()
        elif opt in ("-r", "--rnaInputH5ad"):
            rnaInputH5ad = arg
        elif opt in ("-a", "--atacInputH5ad"):
            atacInputH5ad = arg
        elif opt in ("-w", "--wnnDir"):
            wnnDir = arg
        elif opt in ("-l", "--adjLikelihood"):
            adjLikelihood = True
        elif opt in ("-o", "--outDir"):
            print(arg)
            outDir = arg
  
    adata_rna = sc.read(rnaInputH5ad)
    adata_atac = sc.read(atacInputH5ad)
    scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)
    
    try:
        wnnDir
        nn_idx = np.loadtxt(wnnDir + "/nn_idx.txt", delimiter=',')
        nn_dist = np.loadtxt(wnnDir +"/nn_dist.txt", delimiter=',')
        nn_cells = pd.Index(pd.read_csv(wnnDir +"/nn_cells.txt", header=None)[0])
        mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)
        
    except NameError:
        pass

    adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                        adata_atac,
                                        max_iter=5, 
                                        init_mode="invert", 
                                        verbose=False,
                                        parallel=True,
                                        #n_jobs = 15,
                                        save_plot=False,
                                        rna_only=False,
                                        fit=True,
                                        n_anchors=500 
                                        )
                                        
    if adjLikelihood:
         mv.set_velocity_genes(adata_result, likelihood_lower=0.02)
         
         
    mv.velocity_graph(adata_result)
    mv.latent_time(adata_result)
    
    adata_result.write(outDir+"/multivelo_result.h5ad")

if __name__ == "__main__":
    main(sys.argv[1:])
    






