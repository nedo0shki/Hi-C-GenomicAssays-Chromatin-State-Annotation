# This script gets bedpe file of all contacts and resolution, and
# return COO file of inter-chromosomal contacts and meta data file for mapping
# node IDs and their corresponding position in genome
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
    chr_pos_dict = {}
    for chr_name in valid_chrom:
        num_bin = math.ceil(chrom_size[chr_name]/res)
        for i in range(num_bin):
            chr_pos_dict[(chr_name, i)] = id
            id = id + 1
    return(chr_pos_dict)


def bedpe_to_coo(bedpe_file, res, id_dict, output):

    f = open(output, "w")
    with open(bedpe_file, "r") as bedpe:
        for line in bedpe:
            chr1, start1, end1, chr2, start2, end2, value = line.split()
            first_id = id_dict[(chr1, int(int(start1)/res))]
            second_id = id_dict[(chr2, int(int(start2)/res))]
            f.write(str(first_id) + "\t" + str(second_id) + "\t" + value + "\n")
            if first_id != second_id:
                f.write(str(second_id) + "\t" + str(first_id) + "\t" + value + "\n")
    f.close()


def main():

    parser = argparse.ArgumentParser(description='convert bedpe file to coo file as a LINE input file')
    parser.add_argument('-i', "--bedpe_file_path", required = True,
                       help='Path of Hi-C bedpe file')
    parser.add_argument('-s', '--chrom_size', required = True,
                       help='Path of a file including chromosomes size')
    parser.add_argument('-r', '--resolution', required = True, type = int,
                       help='resolution of annotation')
    parser.add_argument('-o', '--output_file_path', required = True,
                       help='Path of output file to save COO')
    args = parser.parse_args()


    HiC_bedpe_path = args.bedpe_file_path
    chrom_size_file_path = args.chrom_size
    resolution = args.resolution
    output_file = args.output_file_path

    chr_size = read_chrom_size(chrom_size_file_path)
    id_table = make_chr_dict(chr_size, resolution, valid_chrom)
    bedpe_to_coo(HiC_bedpe_path, resolution, id_table, output_file)

if __name__ == "__main__":
    main()
