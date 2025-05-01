import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import scanpy as sc
import os
import sys, getopt
from pathlib import Path

def main(argv):
    matrixFile = ''
    outDir = ''
    try:
        opts, args = getopt.getopt(argv,"hi:g:d:r:o:",["inputH5ad=","gamma=","dims=","reductionKey=","outDir="])
    except getopt.GetoptError:
        print('SEACellsCL.py -i <inputH5ad> -g <gamma> -d <dims> -r <reductionKey> -o <outDir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('SEACellsCL.py -i <inputH5ad> -g <gamma> -d <dims> -r <reductionKey> -o <outDir>')
            sys.exit()
        elif opt in ("-i", "--inputH5ad"):
            inputH5ad = arg
        elif opt in ("-g", "--gamma"):
            gamma = float(arg)
        elif opt in ("-d", "--dims"):
            dimStr= arg
        elif opt in ("-r", "--reductionKey"):
            reductionKey= arg
        elif opt in ("-o", "--outDir"):
            print(arg)
            outDir = arg
          
    print('Output dir is "', outDir)
    print('inputH5ad is "', inputH5ad)
    print('gamma is "', gamma)
    print('dims are "', dimStr)
    print('reductionKey is"', reductionKey)
    ad = sc.read_h5ad(inputH5ad)
    
    # remove attributes added by SeuratDisk that cause SEACells to crash
    del ad.obsp
    del ad.uns
            
    dimStrList = dimStr.split(":")
    build_kernel_on = "X_"+ reductionKey
    ad.obsm[build_kernel_on] = ad.obsm[build_kernel_on][:,range(int(dimStrList[0])-1, int(dimStrList[1]))]
    n_SEACells = int(len(ad)/gamma)
    n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells
    
    model = SEACells.core.SEACells(ad, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)
                  
    model.construct_kernel_matrix()
    # M = model.kernel_matrix
    model.initialize_archetypes()
    #SEACells.plot.plot_initialization(ad, model)
    model.fit(min_iter=10, max_iter=150)
    #model.plot_convergence()
    ad.obs["SEACell"].to_csv(outDir+"/seacellMemberships.csv")

if __name__ == "__main__":
    main(sys.argv[1:])
    






