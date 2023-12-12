##########purpose: take only smallest maf and get gene-level ino
##########input: ind maf files, params of X percent difference of MAF F/M
##########output: 
############(1) file of MAF difference between M/F of more than X percent
############collapsed files for (2) males (3) females and (4) both that is as follows:
############ when referring to mafs, unless explicitly otherwise this refers to the min MAF
############gene/ count of MAF<0.01 for ind /  MAF GTEX/ MAF GNOMAD nonFIN BOTH / 
############MAF GNOMAD nonFIN F/ MAF GNOMAD nonFIN M / minChr / Min Pos /
############ minRef / MinALT / minRefALT / min percent differ / geneTYPE
import argparse
import glob
import re
import gzip

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--infile', type=str, help='directory to read file from')
parser.add_argument('--level', type=str, help='either enst or ensg')
parser.add_argument('--out', type=str, help='outfile collaped ')
args = parser.parse_args()

infile=args.infile
outfile=args.out
level=args.level
print("arguments parsed.")




out=open(outfile,"w")

#keys is gene, value z score
top_gene_dic=dict()
count=0
#index for store line list to use
#if(level=="genes"):
#	myindex=9
#elif(level=="isoforms"):
#	myindex=10
#else:
#	exit("ERROR: INCORRECT INDEX -- NOT GENES OR ISOFORMS")
with open(infile,"r") as f_read:
	for line in f_read.readlines():
		# line_split=(line.decode('utf-8')).split("\t")
		line_split=(line).split("\t")
		#get interesting columns
		count+=1
		chrom=line_split[0]
		start=line_split[1]
		end=line_split[2]
		cluster=line_split[3]
		sample=line_split[4].split('.')[0]
		z=line_split[5]
		# rank_overall=line_split[6] #chr
		# rank_sample=line_split[7] #start
		# rank_cluster=line_split[8] #stop
		
		end_of_line=line_split[-1]
		#print(end_of_line)
		#print(end_of_line.split('gene_id \"')[1].split('"')[0])
		gene_spl=end_of_line.split('gene_id \"')
		transcript="NA"
		if(len(gene_spl)>1):
			gene=(gene_spl[1].split('"')[0])
		#transcript=(end_of_line.split('transcript_id \"')[1].split('"')[0])
		if(level=="genes"):
			enst_or_ensg=gene
		if(level=="isoforms"):
			enst_or_ensg="IN PROGRESS" #transcript
		#(enst_or_ensg)
		store_line=[chrom,start,end,cluster, sample,z,transcript,gene]

		#if there is already a RV recorded for this gene
		if enst_or_ensg not in top_gene_dic:
			top_gene_dic[enst_or_ensg]=store_line
		else:
			stored_z=top_gene_dic[enst_or_ensg][5] #current top z score
			#if z score in this line is more than the best one for this gene/transcript, use this line instead
			if(abs(float(z))>abs(float(stored_z))):
					top_gene_dic[enst_or_ensg]=store_line
	##writes the dictionary to a file
	for this_enst_or_ensg in top_gene_dic:
		this_line=top_gene_dic[this_enst_or_ensg]+['\n']
		out.write(('\t'.join(map(str,this_line))))

out.close()
