import argparse
import glob
import re
import gzip

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--gtf', type=str, help='gtf infile to read from')
parser.add_argument('--outfile', type=str, help='outfile')
args = parser.parse_args()

gtf_infile=args.gtf
outfile=args.outfile
outfile_w=gzip.open(outfile,"wb")

#gtf="/oak/stanford/groups/smontgom/shared/UDN/Tools/.vep/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.101.gtf.gz"

#for speed start with compiling matches you want to make
#some of these have invisible spaces before, so match 0 or 1 spaces in front of this too
re_ensg=re.compile("[ ]{0,1}gene_id")
re_enst=re.compile("[ ]{0,1}transcript_id")
re_hg_name=re.compile("[ ]{0,1}gene_name")
re_genetype=re.compile("[ ]{0,1}gene_biotype")

#in case this category does not exist (which is true for many lines)
#this will return NA instead of breaking with an error
def try_except(re_cat,spl_line):
	found_re=([col for col in spl_line if re_cat.match(col)])
	try:
		re_match=(found_re[0]).split('"')[1]
		return re_match
	except:
		re_match="NA"
		return re_match

		


with gzip.open(gtf_infile,"r") as gtf_read:
	for line in gtf_read.readlines():
		line_split=(line.decode('utf-8')).split("\t")
		chr=line_split[0]
		category=line_split[2] #exon, gene, or transcript		
		start=str(line_split[3])
		end=str(line_split[4])
		details=line_split[8]
		details_split=details.split(';')
		ensg=try_except(re_ensg,details_split)
		enst=try_except(re_enst,details_split)
		hg_name=try_except(re_hg_name,details_split)
		genetype=try_except(re_genetype,details_split)
		outfile_w.write(('\t'.join([chr,start,end,category,ensg,enst,hg_name,genetype,'\n'])).encode())
outfile_w.close()


