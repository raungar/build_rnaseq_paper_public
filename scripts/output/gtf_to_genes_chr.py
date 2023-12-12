import re 
import glob 
import argparse 
import gzip



parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--gtf', type=str, help='gtf file')
parser.add_argument('--infile', type=str, help='gtf file')
parser.add_argument('--outfile', type=str, help='gtf file')
args = parser.parse_args()
gtf=args.gtf
used_genes=args.infile
outfile=args.outfile


#gtf="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf"
#used_genes="Output/sexdeg_v8/all_genes.txt"
#


#print gtf files which are 2 cols: 1=gene_name, 2=[protein_coding/lincRNA] sorted by gene name for aut and x separately
#less -S $gtf | awk -F"[\t;]" '{if($1~/chr[0-9]{1,2}/){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' | sort -k 1 > ${out_a}
#less -S $gtf | awk -F"[\t;]" '{if($1=="chrX"){print $10"\t"$11}}' | sed 's/\"//g' | awk -F" " '{if($4=="protein_coding" || $4=="lincRNA"){print $2"\t"$4}}' | sort -k 1 > ${out_x}

chr_dic={}

with open(gtf,"r") as f_read:
	for line in f_read.readlines():
		if (line.startswith("#")):
			continue
		line_spl=re.split("\t|;",line)
		chr=line_spl[0]
		gene_only=((line_spl[9]).split("\""))[1]
		chr_dic[gene_only]=chr	


outfile_wr=open(outfile,"w")
with open(used_genes,"r") as f_read:
	for l in f_read.readlines():
		gene=l.strip()
		outfile_wr.write(chr_dic[gene]+"\t"+gene+"\n")

outfile_wr.close()
