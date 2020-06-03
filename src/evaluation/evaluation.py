import os
import numpy as np
import statistics as st

def signal_variance_explained(annotation_bedgraph, signal_wig, resolution, annot_vir_res):
    signal_wig = open(signal_wig, "r")
    annotation_bedgraph = open(annotation_bedgraph, "r")
    signals = [[] for i in range(22)]
    last_chr_num = 0
    average = 0
    average_coverage = 0
    for line in signal_wig:
        if not line.startswith("#"):
            chr_name, start, end, value = line.split()
            start = int(start)
            end = int(end)
            value = float(value)
            chr_num = chr_name[3:]
            if str.isdigit(chr_num):
                chr_num = int(chr_num)
                if chr_num != last_chr_num:
                    if last_chr_num != 0:
                        signals[last_chr_num-1].append(average)
                    average = 0
                    average_coverage = 0
                coverage = end - start
                while average_coverage + coverage >= resolution:
                    average = ((average * average_coverage) + (value *
                                                               (resolution - average_coverage))) / resolution
                    signals[chr_num-1].append(average)
                    #virtual_wig.write(str(np.arcsinh(average)) + "\n")
                    average = value
                    coverage = coverage - (resolution - average_coverage)
                    average_coverage = 0
                else:
                    if average_coverage + coverage == 0:
                        average_coverage = 0
                        average = 0
                    else:
                        average = ((average * average_coverage) +
                                   (value * coverage)) / (average_coverage + coverage)
                        average_coverage = average_coverage + coverage
                last_chr_num = chr_num
    signal_dict = {}
    for line in annotation_bedgraph:
        if not line.startswith("#"):
            chr_name, start, end, label = line.split()[0:4]
            start = int(start)
            end = int(end)
            if not annot_vir_res:
                start = int(start/resolution)
                end = int(end/resolution)
            chr_num = chr_name[3:]
            if str.isdigit(chr_num):
                chr_num = int(chr_num)
                if label in signal_dict:
                    signal_dict[label].extend(signals[chr_num-1][start:end])
                else:
                    signal_dict[label] = signals[chr_num-1][start:end]

    signal_mean = {}
    d = []
    x = []
    for key in signal_dict:
        signal_mean = np.mean(signal_dict[key])
        d.extend([s - signal_mean for s in signal_dict[key]])
        x.extend(signal_dict[key])
    print(np.mean(d))
    ve = st.stdev(x) - st.stdev(d)
    print(ve)


signal_variance_explained("annot_test.bedgraph","test.wig",10000,0)
