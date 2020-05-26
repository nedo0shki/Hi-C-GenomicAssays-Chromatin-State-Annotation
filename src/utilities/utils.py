import os
import math
import numpy as np
from scipy.io import loadmat

chrom_sizes = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,
59128983,63025520,48129895,51304566]




# This function create interchromosomal interaction COO from all files in starting with input
# in dump directory. Dump directory could be made using function hicToMat (ref: from SNIPER)

def node_num(chr_num, pos, resolution):
    acc_nums = [math.ceil(x/resolution) for x in chrom_sizes[0:chr_num-1]]
    s = sum(acc_nums)
    s = s + math.ceil(pos/resolution)
    return(s)

def construct_interchrom_COO(dump_dir, resolution):
    o = open("../../data/HiC/interchrom_coo.txt", "a")
    valid_files = [filename for filename in os.listdir(dump_dir) if filename.startswith("input") & filename.endswith(".txt")]
    for file in valid_files:
        file_basename = os.path.splitext(file)[0]
        _, first_chrom_name, second_chrom_name = file_basename.split("_")
        first_chr_num = int(first_chrom_name[4:])
        second_chr_num = int(second_chrom_name[4:])
        file_path = os.path.join(dump_dir, file)
        with open(file_path) as f:
            for line in f:
                first_chr_pos, second_chr_pos, value = line.split()
                first_node_num = node_num(first_chr_num, int(first_chr_pos), resolution)
                second_node_num = node_num(second_chr_num, int(second_chr_pos), resolution)
                o.write(str(first_node_num) + " " + str(second_node_num) + " " + value + "\n")
    o.close()

# ref: from SNIPER
def hicToMat(h,juicer_path,tmp_dir='.',prefix='hic',autoremove=False,overwrite=False,save_matrix=False,verbose=False):
	try:
		os.stat(tmp_dir)
	except:
		os.mkdir(tmp_dir)

	""" Calls juicer_tools to extract hic data into txt files """
	for chrm1 in range(1,23,2):
		for chrm2 in range(2,23,2):
			output_path = os.path.join(tmp_dir,'{2}_chrm{0}_chrm{1}.txt'.format(chrm1,chrm2,prefix))

			if os.path.isfile(output_path) and not overwrite:
				continue

			cmd = 'java -jar {0} dump observed KR {1} {2} {3} BP 100000 {4} > tmp_juicer_log'.format(juicer_path,h,chrm1,chrm2,output_path)
			call([cmd],shell=True)

	""" File path of the inter-chromosomal matrix """
	M_filepath = os.path.join(tmp_dir,'%s_matrix.mat' % prefix)

	""" If overwrite flag is set, reconstruct matrix and check save flag """
	if overwrite:
		M = constructAndSave(tmp_dir,prefix)
	else:
		""" If overwrite unset, check if filepath exists and either load mat or construct new """
		if os.path.isfile(M_filepath):
			M = loadmat(M_filepath)['inter_matrix']
		else:
			M = constructAndSave(tmp_dir,prefix)

	""" If autoremove is set, remove hic .txt files """
	if autoremove:
		for chrm1 in range(1,23,2):
			for chrm2 in range(2,23,2):
				file_path = os.path.join(tmp_dir,'{2}_chrm{0}_chrm{1}.txt'.format(chrm1,chrm2,prefix))
				try:
					os.remove(file_path)
				except:
					continue
	return M

# This function convert node2vec output format (first line number of LVs and number of nodes,
# other lines embedding of a node per line) to wigfix file format
def node2vec_output_to_wigfix(input_file_path, output_dir_path):
    result = []
    with open(input_file_path) as f:
        header_line = next(f)
        num_embed = int(header_line.split()[1])
        for line in f:
            line = line.split()
            pos = line[0] #This will get the first 3 characters of the line
            pos = int(pos)
            thisLine = {"pos":pos, "embedding":line[1:]}
            result.append(thisLine)
    sortedList = sorted(result, key=lambda k: k["pos"])
    print(sortedList[0:10])
    chr_num = 21
    resolution = 1
    i = 0
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path)
    embed = []
    for e in range(num_embed):
        embed.append(open(os.path.join(output_dir_path, "embedding" + str(e+1) + ".wigfix"), "w"))
    while i < len(sortedList):
        header = ("fixedStep chrom=chr" + str(chr_num) + " start=" +
         str(resolution*sortedList[i]['pos']) + " step=" + str(resolution) +
          " span=" + str(resolution))
        for e in range(num_embed):
            embed[e].write(header)
            embed[e].write("\n")
            embed[e].write(sortedList[i]['embedding'][e])
            embed[e].write("\n")
        i = i + 1
        while sortedList[i]['pos'] == sortedList[i-1]['pos']+1:
            for e in range(num_embed):
                embed[e].write(sortedList[i]['embedding'][e])
                embed[e].write("\n")
            i = i + 1
            if i == len(sortedList):
                break

# This function gets a wig file and resolution and output new wigfix file each bin representing
# average value over $resolution bins.
def make_vir_res(wig_file_path, resolution):
    resolution = int(resolution)
    original_wig = open(wig_file_path, "r")
    wig_path, wig_name = os.path.split(wig_file_path)
    wig_name, _ = os.path.splitext(wig_name)
    wig_name = wig_name + "-VirRes" + str(resolution) + ".wigfix"
    vir_wig_path = os.path.join(wig_path, wig_name)
    virtual_wig = open(vir_wig_path, "w")
    last_chr_name = ""
    average = 0
    average_coverage = 0
    for line in original_wig:
        if not line.startswith("#"):
            chr_name, start, end, value = line.split()
            start = int(start)
            end = int(end)
            value = float(value)
            if chr_name != last_chr_name:
                if last_chr_name != "":
                    virtual_wig.write(str(average) + "\n")
                average = 0
                average_coverage = 0
                virtual_wig.write("fixedStep chrom=" + chr_name + " start=1 step=1 span=1" + "\n")
            coverage = end - start
            while average_coverage + coverage >= resolution:
                average = ((average * average_coverage) + (value * (resolution - average_coverage))) / resolution
                virtual_wig.write(str(np.arcsinh(average)) + "\n")
                average = value
                coverage = coverage - (resolution - average_coverage)
                average_coverage = 0
            else:
                if average_coverage + coverage == 0:
                    average_coverage = 0
                    average = 0
                else:
                    average = ((average * average_coverage) + (value * coverage)) / (average_coverage + coverage)
                    average_coverage = average_coverage + coverage
            last_chr_name = chr_name

# This function gets outputs of SNIPER embedding (embedding mat file and index file) and the
# resolution of embeddings (if we want to have virtual resolution equal to resolution of SNIPER,
# we should set this parameter to 1) and save a wigfix file per latent variable in a wigfix_dir_path

def SNIPER_output_to_wigfix(data_path, index_path, wigfix_dir_path, resolution):
    data_mat = loadmat(data_path)
    index_mat = loadmat(index_path)
    num_LV = data_mat['odd_encoded'].shape[1]
    if not os.path.exists(wigfix_dir_path):
        os.makedirs(wigfix_dir_path)
    LV = []
    for l in range(num_LV):
        LV.append(open(os.path.join(wigfix_dir_path, "LV" + str(l+1) + ".wigfix"), "w"))
    i = 0
    while i < index_mat['rowMap'].shape[0]:
        header = "fixedStep chrom=chr" + str(index_mat['rowMap'][i,1]) + " start=" + str((resolution*index_mat['rowMap'][i,2])+1) + " step=" + str(resolution) + " span=" + str(resolution)
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
        header = "fixedStep chrom=chr" + str(index_mat['colMap'][i,1]) + " start=" + str((resolution*index_mat['colMap'][i,2])+1) + " step=" + str(resolution) + " span=" + str(resolution)
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
