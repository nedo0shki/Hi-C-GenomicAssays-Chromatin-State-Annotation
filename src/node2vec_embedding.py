# input is a hic file and output is a directory including an embedding file per chromosome

import os
import sys
import pandas
import numpy
import argparse
import utils





def main():
    parser = argparse.ArgumentParser(description='Make dir of correlation matrix from .hic file')
    parser.add_argument('-h', "--hic_file", required = True,
                       help='hic file path')
    parser.add_argument('-r', '--resolution', required = True, type = int,
                       help='resolution of matrices')
    parser.add_argument('-d', '--dump_dir', required = True,
                       help='path of dump directory to store matrices in')
    parser.add_argument('-j', '--juicer', required = True,
                       help='path of juicer tool')
    args = parser.parse_args()

    hic_path = args['hic_file']
    resolution = args['resolution']
    dump_dir = args['dump_dir']
    juicer_path = args['juicer']

    utils.hicToCorrIntraMat(hic_path, dump_dir, resolution, juicer_path)
    


if __name__ == "__main__":
    main()
