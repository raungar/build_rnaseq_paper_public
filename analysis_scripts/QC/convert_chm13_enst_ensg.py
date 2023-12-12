import re

myfile="/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/CHM13v2.0/chm13.draft_v2.0.gene_annotation.gff3"
outfile="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/chm13_enstANDensg.txt"
my_out=open(outfile,"w")
with open(myfile) as f:
    for line in f:
        l_spl=line.split("\t")
        if(len(l_spl)==1):
            continue
        if(l_spl[2]=="transcript"):
            long_text=l_spl[8]
            enst=re.search("source_transcript=([^;]*);",long_text)
            ensg=re.search("source_gene=([^;]*);",long_text)
            gene_chm13=re.search("gene_id=([^;]*);",long_text)
            transcript_chm13=re.search("transcript_id=([^;]*);",long_text)
            my_line=[ensg.group(1),enst.group(1),gene_chm13.group(1),transcript_chm13.group(1),'\n']
            my_out.write('\t'.join(my_line))
