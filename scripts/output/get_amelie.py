  
import sys, gzip, csv, numpy, os, requests, urllib3, subprocess, ast, argparse, operator
from itertools import chain

def get_args():
	parser = argparse.ArgumentParser(prog='VCF Annotator', usage='%(prog)s [options]')
	parser.add_argument("--file",type=str,help="Proband File Path",required=True)
	parser.add_argument("--metadata",type=str,help="metadata fiel",required=True)
	parser.add_argument("--sample",type=str,help="sample id",required=True)
	parser.add_argument("--out",type=str,help="Output File Name",required=True)
	args=parser.parse_args()
	return(args)

#query amelie for genes in the line against the phenotype and get score
#returns the gene(s) queried, the list of scores, and the top score
def get_amelie(input_genes,phenotypes):

	#submit request
	r = requests.post('https://medgen.stanford.edu/api/',
		verify=False,
		auth = ("udn", "udnAmelie"),
		data = {
			"genes" : input_genes,
			"phenotypes" : phenotypes
		}
	)

	amelie_request=ast.literal_eval((r.text).strip())

	#print("AMELIE genes then phenotypes")
	#print(input_genes)
	#print(phenotypes)
	#print("REQUEST")
	#print(amelie_request)

	#Set to none in case nothing found
	amelie_gene="None"
	amelie_max="None"
	amelie_list="None"

	#If request found, parse through this list of lists
	if amelie_request:
		amelie_gene=set()
		for sublist in amelie_request:
			#check if this is the gene name or a list of scores
			for item in sublist:
				#if is a list
				if isinstance(item,list):
					#check if score exists
					if not item:
						amelie_max="None"
						amelie_list=["None"]
					else:
						#find the max score, and add all scores to the score list
						amelie_max=0
						amelie_list=[]
						for element in item:
							amelie_list.append(element[0])
							if(float(element[0])>amelie_max):
								amelie_max=float(element[0])
				else:
					#then this is the gene, so add it to the list
					amelie_gene.add(item.strip())
		#reformat for return
		amelie_gene=';'.join(amelie_gene).strip()
		amelie_list=';'.join(map(str,amelie_list)).strip()
	return amelie_gene, amelie_max, amelie_list


def main():
	args=get_args()
		
	##open output
	##read in file	
		#get hpo terms
		#read in line
		#read in gene name
		amelie_gene, amelie_max, amelie_list=get_amelie(line_gene_names,hpo)

if __name__ == "__main__":
	main()
