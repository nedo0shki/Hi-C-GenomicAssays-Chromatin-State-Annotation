# This script gets embedding file with node id and chromosome size file, and return
# embeding file including chromosome name and position corresponding to ids

import math
import argparse

valid_chrom = ["chr" + str(i+1) for i in range(22)]

def read_chrom_size(infile):

    chrom_size = {}
    with open(infile, "r") as chr_size:
        for line in chr_size:
            name, size = line.split()
            chrom_size[name] = int(size)
    return(chrom_size)

def make_chr_dict(chrom_size, res, valid_chrom):

    id = 0
    chr_name_dict = {}
    chr_pos_dict = {}
    for chr_name in valid_chrom:
        num_bin = math.ceil(chrom_size[chr_name]/res)
        for i in range(num_bin):
            chr_pos_dict[id] = i
            chr_name_dict[id] = chr_name
            id = id + 1
    return([chr_name_dict, chr_pos_dict])

def set_index(embedding_file, res, chr_name_dict, chr_pos_dict, output):

    f = open(output, "w")
    with open(embedding_file, "r") as embedding:
        for line in embedding:
            l = line.split()
            id = int(l[0])
            emb = l[1:]
            chr_name = chr_name_dict[id]
            chr_pos = chr_pos_dict[id]
            f.write(chr_name + " " + str(chr_pos*res) + " " + str((chr_pos+1)*res))
            for e in emb:
                f.write(" " + e)
            f.write("\n")
    f.close()

def main():

    parser = argparse.ArgumentParser(description='lift bedpe file to other assembly')
    parser.add_argument('-i', "--embedding_file_path", required = True,
                       help='Path of embedding file')
    parser.add_argument('-s', '--chrom_size', required = True,
                       help='Path of a file including chromosomes size')
    parser.add_argument('-r', '--resolution', required = True, type = int,
                       help='resolution of annotation')
    parser.add_argument('-o', '--output_file_path', required = True,
                       help='Path of output file to save indexed embedding')
    args = parser.parse_args()


    embedding_path = args.embedding_file_path
    chrom_size_file_path = args.chrom_size
    resolution = args.resolution
    output_file = args.output_file_path

    chr_size = read_chrom_size(chrom_size_file_path)
    chr_name_table, chr_pos_table = make_chr_dict(chr_size, resolution, valid_chrom)
    set_index(embedding_path, resolution, chr_name_table, chr_pos_table, output_file)

if __name__ == "__main__":
    main()
