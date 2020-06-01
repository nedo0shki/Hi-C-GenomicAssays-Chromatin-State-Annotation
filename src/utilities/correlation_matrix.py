import os
import sys
import pandas
import numpy
import argparse

parser = argparse.ArgumentParser(description='Make dir of correlation matrix from .hic file')
parser.add_argument('-h', "--hic_file", required = True,
                   help='hic file path')
parser.add_argument('-r', '--resolution', required = True, type = int,
                   help='resolution of matrices')
parser.add_argument('-d', '--dump_dir', required = True,
                   help='path of dump directory to store matrices in')
args = parser.parse_args()

hic_path = args['hic_file']
resolution = args['resolution']
dump_dir = args['dump_dir']
