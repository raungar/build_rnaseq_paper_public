import argparse
import numpy as np
import random
from sklearn import preprocessing
import pandas as pd
from timeit import default_timer as timer
from itertools import chain


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile',help='input file', required=True)
    #infile: /oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt
    parser.add_argument('--outfile',help='output file',required=True)
    parser.add_argument('--alpha',help='alpha value for trimmed mean, default is 0.1',default=0.05)
    parser.add_argument('--num_bootstrap',help='number of B bootstraps',required=True)
    args = parser.parse_args()
    return(args)

def get_trimmed_mean(samples, alpha):
    sample_len=len(samples)
    num_samples_to_drop=np.floor(alpha*sample_len).astype(int)
    samples_sorted=np.sort(samples)
    samples_trimmed=samples_sorted[num_samples_to_drop:sample_len-num_samples_to_drop]
    return(np.mean(samples_trimmed))

def bootstrap(xs,num_b,alpha):
    trimmed_means=[]

    for i in range(int(num_b)):
        randomly_sampled_xs=random.choices(xs,k=len(xs)) #randomly sample with replacement
        trimmed_means.append(get_trimmed_mean(randomly_sampled_xs,alpha))
    return(trimmed_means)

def get_new_zs(bootstrapped_gene_means,xs):
    bootstrapped_gene_std=np.std(bootstrapped_gene_means)
    bootstrapped_gene_means_mean=np.mean(bootstrapped_gene_means)
    new_zs=(xs-bootstrapped_gene_means_mean)/bootstrapped_gene_std
    #print(bootstrapped_gene_means_mean)
    #print(bootstrapped_gene_std)
    return(new_zs)


def main():
    args=parse_args()
    print(args.infile)

    df=pd.read_csv(args.infile,  sep='\t')
    num_uniq_genes=len(df.gene.unique())
    bootstrapped_list=[]
    for this_gene in df.gene.unique():
        start=timer()
        print(this_gene)
        this_gene_df=df.loc[df['gene']==this_gene]
        xs=this_gene_df.zscore.tolist()
        bootstrapped_gene_means=bootstrap(xs,args.num_bootstrap,args.alpha)
        new_zs=get_new_zs(bootstrapped_gene_means,xs)
        new_zs_scaled=preprocessing.scale(new_zs)
        print(xs)
        print(new_zs)
        print(pd.DataFrame(list(zip(xs,new_zs_scaled,new_zs))))
        bootstrapped_list.append(new_zs) #unlist

        end=timer()
        print("time: "+str(end-start)+". This is x "+str(num_uniq_genes)+ " genes so it'll take "+str(num_uniq_genes*(end-start))+" secs")
        break
    print(list(chain.from_iterable(bootstrapped_list)))


main()

