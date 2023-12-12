#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import gzip
import numpy as np
import os
############## FUNCTIONS

def sample_dict_maker(sample_fn):
	'''
	Reads in sample annotation file and makes dict. Keys are sample IDs, and values are tissues.
	'''
	sample_fh = open(sample_fn, 'r')
	sample_dict = {}
	for line in sample_fh:
		sample, tissue = line.rstrip().split('\t')
		if tissue != '':
			sample_dict[sample] = tissue
	sample_fh.close()
	return sample_dict

def split_by_tissue(sample_dict, gtex_fn, outfile,this_tissue):
	'''
	Read in GTEx TPM or reads file and split by tissues.
	'''
	gtex = gzip.open(gtex_fn, 'r')
	for line in gtex:
		line = line.decode('utf8').rstrip().split('\t')
		if 'Name' in line:
			header=line
			break
		if 'gene_id' in line:
			header = line
			break

	tissue_dict = {}
	for sample in header[2:]:
		if sample in sample_dict:
			tissue = sample_dict[sample]
			if tissue not in tissue_dict:
				tissue_dict[tissue] = []
			tissue_dict[tissue].append(header.index(sample))
	else:
		print(sample, 'not in sample-to-tissue map. Skipping it.', file = sys.stderr)

	#print("this is tissue dict")
	#print(tissue_dict)
	header = np.array(header)
	file_dict = {}
	for tissue in tissue_dict:
		if tissue != this_tissue:
			continue
		print("for tissue in tissue dict")
		file_dict[tissue] = open(outfile, 'w')
		print("open")
		#print(header)
		out_names = header[tissue_dict[tissue]]
		for i in range(len(out_names)):
			out_names[i] = '-'.join(out_names[i].split('-')[0:2])
		print('\t'.join(['Gene'] + out_names.tolist()), file = file_dict[tissue])
	print("for line")
	for line in gtex:
		line = np.array(line.decode('utf8').rstrip().split())
		for tissue in tissue_dict:
			if tissue != this_tissue:
				continue
			print('\t'.join([line[0]] + line[tissue_dict[tissue]].tolist()), file = file_dict[tissue])
	print("for tiss in file dict")
	for tissue in file_dict:
		if tissue != this_tissue:
			continue

		file_dict[tissue].close()



############## MAIN
print("split_expr_by_tissues.py")

usage = 'Splits combined RPKM file from GTEx into expression matrices by tissue.'

parser = argparse.ArgumentParser(description = usage)
parser.add_argument('--gtex', required = True, help = 'GTEx combined TPM or read count file')
parser.add_argument('--outfile', required = True, help = 'outfile path')
parser.add_argument('--sample', required = True, help = 'sample annotation file')
parser.add_argument('--tissue', required = True, help = 'tissue to grab')
##parser.add_argument('--end', required = True, help = 'file ending')

args = parser.parse_args()

###RAU at out if does not exist
#if not os.path.exists(args.out):
#	os.makedirs(args.out)

## Make sample dict and get list of tissues
sample_dict = sample_dict_maker(args.sample)

print(sample_dict)

## Split GTEx combined file by tissue
split_by_tissue(sample_dict, args.gtex, args.outfile,args.tissue)
