# this script gets an embedding file and generate some analysis results

import pandas as pd
import numpy as np
import sklearn
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import math
import statistics as st
import argparse
import matplotlib.backends.backend_pdf
import os



parser = argparse.ArgumentParser(description='Generate signals enrichment heatmap')
parser.add_argument('-i', "--input_file", required = True,
                   help='embedding file path')
parser.add_argument('-o', '--output', required = True,
                   help='Path of Pdf to save results in')
args = parser.parse_args()


def read_annot_file(annot_path, vir_res, resolution):
    annot_df = pd.read_csv(annot_path, sep = "\t", header = None)
    annot_df = annot_df.iloc[:,0:4]
    if not vir_res:
        annot_df.iloc[:,1:3] = (annot_df.iloc[:,1:3]/resolution).astype(int)
    annot_df.columns = ['chr_name', 'start', 'end', 'label']
    return(annot_df)

def melt_annotation(annot_df):
    annot_df = annot_df.reset_index()
    melt = annot_df.melt(id_vars=['index', 'chr_name', 'label'], value_name = 'pos').drop('variable', axis=1)
    melt_annot_df = melt.groupby('index').apply(lambda x: x.set_index('pos', drop = False)\
                                    .reindex(range(x.loc[:,'pos'].values[0], x.loc[:,'pos'].values[1])))\
                                .ffill().drop(['pos','index'], axis = 1).reset_index(level = 'pos').reset_index(drop=True)
    return(melt_annot_df)

def read_embedding_file(embedding_path, scaled):

    emb = pd.read_csv(embedding_path, sep = " ", header = None)
    emb.iloc[:,1] = (emb.iloc[:,1]/100000).astype(int)
    emb = emb.sort_values(by = [0,1])
    if scaled:
        emb.iloc[:,3:len(emb.columns)] = StandardScaler().fit_transform(emb.iloc[:,3:len(emb.columns)])
    for i in range(3,len(emb.columns)):
        emb = emb.rename(columns = {i: "emb" + str(i-2)})
    features = emb.iloc[:,3:len(emb.columns)]
    kmeans = KMeans(n_clusters=5)
    kmeans.fit(features)
    emb['kmeans_label'] = kmeans.predict(features)
    reducer = umap.UMAP()
    umaps = reducer.fit_transform(features)
    emb['umap1'] = umaps[:,0]
    emb['umap2'] = umaps[:,1]
    emb.drop(emb.columns[2], axis=1, inplace=True)
    emb = emb.rename(columns = {0:'chr_name', 1:'pos'})
    return(emb)

def variance_explained(signal_vector, label_vector):

    d = []
    x = []
    for l in np.unique(label_vector):
        curr_x = signal_vector[label_vector==l]
        pred = np.mean(curr_x)
        curr_d = [cx - pred for cx in curr_x]
        d.extend(curr_d)
        x.extend(curr_x)
    ve = (st.stdev(x) - st.stdev(d))/st.stdev(x)
    return(ve)

def set_int_value(labels):
    value_to_int = {value: i for i, value in enumerate(sorted(pd.unique(labels)))}
    int_labels = [value_to_int[i] for i in labels]
    int_labels = np.asarray(int_labels).reshape(1,-1)
    return(int_labels)

def make_signals_df(signals_folder_path):
    signals_df = pd.DataFrame()
    for file in os.listdir(signals_folder_path):
        signal = pd.read_csv(os.path.join(signals_folder_path, file), sep = "\t", header = None)
        signal = signal.iloc[:,[0,1,3]]
        signal_name = os.path.splitext(file)[0]
        signal.columns = ["chr_name", "pos", signal_name]
        if signals_df.empty:
            signals_df = signal
        else:
            signals_df = pd.merge(signals_df, signal, how = 'inner', on = ['chr_name','pos'])
    return(signals_df)

def main():
    pdf = matplotlib.backends.backend_pdf.PdfPages(args.output)
    SC_annotation = read_annot_file("../../../data/annotations/Rao_SC.bed", 1, 100000)
    SC_annotation = melt_annotation(SC_annotation)
    SC_annotation.rename(columns = {'label': 'SC'}, inplace = True)
    RT = pd.read_csv("../../../data/RT/six_phase/RT_res100000.txt", sep = "\t")
    TSA_seq = pd.read_csv("../../../data/genomic-assays/TSA-seq/LaminA_TSA-seq_binned100000.bedgraph", sep = "\t", header = None)
    TSA_seq = TSA_seq.drop(2,axis=1)
    TSA_seq.iloc[:,1] = (TSA_seq.iloc[:,1]/100000).astype(int)
    TSA_seq.columns = ["chr_name", "pos", "Lamin_TSA_seq"]
    embedding = read_embedding_file(args.input_file, 0)
    emb_size = len(embedding.columns)-5
    embedding = pd.merge(embedding, SC_annotation, how = 'inner', on = ['chr_name','pos'])
    embedding = pd.merge(embedding, RT, how = 'inner', on = ['chr_name','pos'])
    embedding = pd.merge(embedding, TSA_seq, how = 'inner', on = ['chr_name','pos'])
    emb_names = ["emb" + str(i+1) for i in range(emb_size)]

    ############# correlation between embeddings ###########
    fig = plt.figure()
    sns.heatmap(np.corrcoef(embedding.loc[:,emb_names].transpose()), cmap = "Blues", annot = True)
    plt.title("Embeddings correlation")
    pdf.savefig(fig)
    ############# correlation of each embedding with TSA_seq ########
    TSA_seq_corr = [pearsonr(embedding['Lamin_TSA_seq'], embedding[emb_name])[0] for emb_name in emb_names]
    fig = plt.figure()
    plt.bar(emb_names,TSA_seq_corr)
    plt.title("Correlation of embeddings with TSA-seq")
    pdf.savefig(fig)
    ############ distribution of TSA_seq for different labels #########
    fig = plt.figure(figsize=(9,5))
    TSA_seq_VE = variance_explained(embedding['Lamin_TSA_seq'], embedding['kmeans_label'])
    for i, group in embedding.groupby("kmeans_label"):
        sns.distplot(group['Lamin_TSA_seq'],
                     label = str(i) + " (" + str(round(np.mean(group['Lamin_TSA_seq']),2)) + ")")
    plt.legend()
    plt.title("Distribution of LaminA TSA seq signal for \n labels from LINE + kmeans (VE: " + str(round(TSA_seq_VE, 2)) + ")")
    pdf.savefig(fig)
    ############ scatter plot of embedding umap colored by SC labels #############
    fig = plt.figure()
    sns.scatterplot(x="emb1", y="emb2", hue="SC", data=embedding)
    pdf.savefig(fig)
    ############ overlap with SC labels #############
    SC_SI = sklearn.metrics.silhouette_samples(embedding[embedding['SC'].isin(["A1","A2","B1","B2","B3","B4"])].loc[:,emb_names],embedding[embedding['SC'].isin(["A1","A2","B1","B2","B3","B4"])]['SC']).mean()
    a = embedding.groupby(['SC','kmeans_label']).size()
    a_table = a.unstack(level=0)
    fig = plt.figure()
    plt.table(cellText=a_table.values, colLabels=a_table.columns, rowLabels=a_table.index,loc='center')
    plt.axis("off")
    plt.title("Overlap of kmeans annotation and SC (SI mean for embeddings based on SC label: " + str(round(SC_SI,2)) + ")")
    pdf.savefig(fig)
    ############# enrichment plot ################
    signals_df = make_signals_df("../../../data/genomic-assays/ChIP-seq/binned_bedgraphs_vir_res")
    signals_name = signals_df.columns[2:]
    embedding = pd.merge(embedding, signals_df, how = 'inner', on = ['chr_name', 'pos'])
    en = embedding.groupby("kmeans_label")[signals_name].mean()
    mean_en = embedding.loc[:,signals_name].mean(axis=0)
    en = round(en.divide(mean_en, axis=1),2)
    fig = plt.figure(figsize=(8,6))
    sns.heatmap(en, annot = True, cmap="PiYG", center = 1, linewidths=.5, square = True, cbar = False)
    plt.title("Enrichment of genomics assays based on \n LINE+kmeans")
    pdf.savefig(fig)
    ############### coverage plot ###########
    fig = plt.figure()
    sns.barplot(x = np.unique(embedding['kmeans_label'], return_counts = True)[0], y = np.unique(embedding['kmeans_label'], return_counts = True)[1])
    plt.title("LINE+kmeans labels coverage")
    pdf.savefig(fig)
    ############### RT ##############
    phases = ["G1", "S1", "S2", "S3", "S4", "G2"]
    phases_ve = []
    for phase in phases:
            phases_ve.append(variance_explained(embedding[phase], embedding['kmeans_label']))
    fig = plt.figure()
    plt.bar(phases, phases_ve)
    plt.title("Variance explained for different phases based on LINE+kmeans labels")
    pdf.savefig(fig)
    #######chr16-17 corr plot##########
    chrom_size = {"chr16": 90354753, "chr17": 81195210}
    res = 100000
    id = 0
    chr_pos_dict = {}
    for key in chrom_size.keys():
        num_bin = math.ceil(chrom_size[key]/res)
        for i in range(num_bin):
            chr_pos_dict[(key, i)] = id
            id = id + 1
    chr16_17_embedding = embedding[embedding['chr_name'].isin(["chr16", "chr17"])]
    chr16_17_embedding['id'] = [chr_pos_dict[i[0], i[1]]
                             for i in zip(chr16_17_embedding['chr_name'],chr16_17_embedding['pos'])]
    chr16_17_array = np.zeros([1716,emb_size])
    chr16_17_array[chr16_17_embedding['id'],:] = chr16_17_embedding.loc[:,emb_names]
    chr16_17_corr_mat = np.corrcoef(chr16_17_array)
    chr16_17_labels = np.empty(id, dtype='U2')
    chr16_17_labels[chr16_17_embedding['id']] = chr16_17_embedding['kmeans_label'].values
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((11,10), (0,0), colspan=10, rowspan=10)
    ax2 = plt.subplot2grid((11,10), (10,0), colspan=10, rowspan=1)
    sns.heatmap(chr16_17_corr_mat, cmap = "Blues", ax = ax1, cbar = False, xticklabels=False, yticklabels=False)
    sns.heatmap(set_int_value(chr16_17_labels), ax = ax2, cbar = False, xticklabels=False, yticklabels=False)
    pdf.savefig(fig)
    pdf.close()

if __name__ == "__main__":
    main()
