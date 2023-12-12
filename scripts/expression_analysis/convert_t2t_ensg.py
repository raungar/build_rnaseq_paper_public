import re 
import glob 
import argparse 
import gzip



parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--gff', type=str, help='gtf file')
parser.add_argument('--infile', type=str, help='gtf file')
parser.add_argument('--outfile', type=str, help='gtf file')
args = parser.parse_args()
gff=args.gff
infile=args.infile
outfile=args.outfile

gene_dic={} #key t2t id, value is ensembl id
transcript_dic={} # key is t2t id, value is ensembl id

#read gff file into a dictionary to access later
with open(gff,"r") as f_read:
	for line in f_read.readlines():
		if (line.startswith("#")):
			continue
		line_spl=re.split("\t|;",line)
		feature_type=line_spl[2] #gene/transcript/exon/intron/etc
		if(feature_type=="gene"):
			ensg=(line.split("source_gene=")[1]).split(";")[0]
			chm13=(line.split("ID=")[1]).split(";")[0]
			gene_dic[chm13]=ensg
		if(feature_type=="transcript"):
			enst=(line.split("source_transcript=")[1]).split(";")[0]
			chm13_transcript=(line.split("transcript_id=")[1]).split(";")[0]
			transcript_dic[chm13_transcript]=enst

outfile_wr=open(outfile,"w")
#replace t2t id with ensembl id
with open(infile,"r") as f_read:
	header=f_read.readline()
	outfile_wr.write(header)
	for l in f_read.readlines():
		line_spl=l.split("\t")
		ensg=gene_dic[line_spl[0]] #replace in place ensg
		enst_list=[] #transcripts are in a comma separated list so have to loop through
		for transcript in line_spl[1].split(","):
			enst_list.append(transcript_dic[transcript]) #append dictionaried transcripts to list
		enst=','.join(enst_list) #rejoin transcripts with a comma as before and replace in place
		write_line=[ensg,enst]+line_spl[2:]
		outfile_wr.write('\t'.join(write_line)) #write line to new outfile

outfile_wr.close()
