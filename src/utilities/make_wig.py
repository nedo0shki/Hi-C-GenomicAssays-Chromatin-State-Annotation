# this script get output of SNIPER (matrix of latent variables, and index matrix) and return N wigfix files, each of them corresponding with one LV track

import os
from scipy.io import loadmat
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--data_path", help = "Path of encoded data from SNIPER", required = True)
ap.add_argument("--index_path", help = "Path of index data from SNIPER", required = True)
ap.add_argument("--LV_wig_dir_path", help = "Path of directory to save wigfix files of LV tracks", required = True)
ap.add_argument("--resolution", help = "resolution of annotation", required = True, type = int)
args = vars(ap.parse_args())

data_mat = loadmat(args['data_path'])
index_mat = loadmat(args['index_path'])
resolution = args['resolution']
num_LV = data_mat['odd_encoded'].shape[1]
if not os.path.exists(args['LV_wig_dir_path']):
    os.makedirs(args['LV_wig_dir_path'])
LV = []
for l in range(num_LV):
    LV.append(open(os.path.join(args['LV_wig_dir_path'], "LV" + str(l+1) + ".wigfix"), "w"))
i = 0
while i < index_mat['rowMap'].shape[0]:
    header = "fixedStep chrom=chr" + str(index_mat['rowMap'][i,1]) + " start=" + str(resolution*index_mat['rowMap'][i,2]) + " step=" + str(resolution) + " span=" + str(resolution)
    for l in range(num_LV):
        LV[l].write(header)
        LV[l].write("\n")
        LV[l].write(str(data_mat['odd_encoded'][i,l]))
        LV[l].write("\n")
    i = i + 1
    while index_mat['rowMap'][i,2] == index_mat['rowMap'][i-1,2]+1:
        for l in range(num_LV):
            LV[l].write(str(data_mat['odd_encoded'][i,l]))
            LV[l].write("\n")
        i = i + 1
        if i == index_mat['rowMap'].shape[0]:
            break
i = 0
while i < index_mat['colMap'].shape[0]:
    header = "fixedStep chrom=chr" + str(index_mat['colMap'][i,1]) + " start=" + str(resolution*index_mat['colMap'][i,2]) + " step=" + str(resolution) + " span=" + str(resolution)
    for l in range(num_LV):
        LV[l].write(header)
        LV[l].write("\n")
        LV[l].write(str(data_mat['even_encoded'][i,l]))
        LV[l].write("\n")
    i = i + 1
    while index_mat['colMap'][i,2] == index_mat['colMap'][i-1,2]+1:
        for l in range(num_LV):
            LV[l].write(str(data_mat['even_encoded'][i,l]))
            LV[l].write("\n")
        i = i + 1
        if i == index_mat['colMap'].shape[0]:
            break
for l in range(num_LV):
    LV[l].close()
