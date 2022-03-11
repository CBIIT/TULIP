#!/usr/bin/python
## CanClass: A tools for performing Cancer type classification
## This script downloads the model from MoDaC, corresponding to teh cancer type 17 or 32 selection###############
## Check the matrix format for the genes 19k genes######
## Print out the samples names to a file ####
## Transform the matrix and prepare it removing the sample name column [first column] and gene name row [first row] and reset the index ####
## Convert the FPKM-UQ values to TPM value and normalize it ##
## Predict script ###

#### Author: Satish RG #######
#### Author: Sara Jones #######


from __future__ import print_function
import os, warnings 
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from keras.models import load_model
from keras.layers import Input, Dense, Dropout, Activation, Conv1D, MaxPooling1D, Flatten
from keras.models import Sequential, Model, model_from_json, model_from_yaml


def get_arg():
    parse = argparse.ArgumentParser()
    reqargs = parse.add_argument_group('Required')
    reqwdef = parse.add_argument_group('Required with default')
    reqargs.add_argument('-i', '--inputmatrix', type=str, required=True, help='Provide the full path of the expression matrix file (FPKM-UQ) in the required format.')
    reqargs.add_argument('-o', '--outputdir', type=str, required=True, help='Provide the full path to the output directory.')
    reqwdef.add_argument('-t', '--cancertypes', type=str, default = '17', help='provide the # of cancer types 17 or 32, you want to use to classify your data')
    allarg = parse.parse_args()
    return allarg

def matrixprep(matfile, outdir):
    ######################################
    #### Output files ######
    ######################################
    sampname=outdir+"sample_name.txt"
    formatrix = outdir+"inputmatrix.tsv"
    ######################################
    ##### Matrix processing ######
    ######################################
    df_FPKM_UQ = pd.read_csv(matfile, low_memory=False, sep="\t")
    df_FPKM_UQ.shape
    ### submitters_id mapping to project_id
    cols = df_FPKM_UQ.columns[1:].values.tolist()
    type(cols)
    with open(sampname, 'w') as f:
        for item in cols:
            f.write("%s\n" % item) 
    ###Transpose the matrix
    dft_FPKM_UQ = df_FPKM_UQ.T
    dft_FPKM_UQ.shape
    # remove the two two rows and save the output
    dftm_FPKM_UQ = dft_FPKM_UQ.drop(dft_FPKM_UQ.index[0:1], axis=0)
    dftm_FPKM_UQ.shape
    dftm_FPKM_UQ = dftm_FPKM_UQ.reset_index(drop=True)
    # Drop the index (sample name row)
    features = dftm_FPKM_UQ.drop(dftm_FPKM_UQ.index[:0], axis=1)
    # TPM Calculations
    sfeatures = features.div(features.sum(axis=1), axis=0)
    sfeatures = sfeatures * 1000000
    sfeatures1 = sfeatures.astype(np.float64).apply(np.log10)
    sfeatures1[sfeatures1 < 0] = 0
    ###### Writing the formated file ready for CanClass #####
    sfeatures.to_csv(formatrix, sep='\t', index=False)


    

def main():
    arg = get_arg()
    print(arg.inputmatrix+"\t"+arg.outputdir+"\n")
    matrixprep(arg.inputmatrix,arg.outputdir)

if __name__ == "__main__":
    main()