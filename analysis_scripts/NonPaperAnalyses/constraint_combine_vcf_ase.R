library("data.table")
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--vcf"), type="character", default=NULL, help="path for vcf", metavar="character"),
  make_option(c("--status"), type="character", default=NULL, help="case or control", metavar="character"),
  make_option(c("--tiss"), type="character", default=NULL, help="tissue", metavar="character"),
  make_option(c("--spl"), type="character", default=NULL, help="path for leafcuttermd", metavar="character"),
  make_option(c("--out"), type="character", default=NULL, help="path for outfile", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
vcf_f=opt$vcf
spl_f=opt$spl
out=opt$out
status=opt$status
tiss=opt$tiss

vcf=fread(vcf_f)
print(head(vcf))
#colnames(vcf)<-c("chr" ,"pos","pos2","ref","alt", "ind","gene","ensg","enst","maf", "cadd_phred","impactful","impact_score","vartype","dna_change","aa_change","disease","gene_family","num_rv","ensg_withperiod","hgnc","regional_score" )
#print(c("chr" ,"pos","ref","alt", "ind","gene","ensg","enst","maf", "cadd_phred","impactful","impact_score","vartype","dna_change","aa_change","disease","gene_family","regional_score","num_rv" ))
#this_chr,pos,ref,alt,ind,gene,ensg,enst,maf_use,cadd_phred,impactful,highest_impact_score,dna_change,aa_change,disease,gene_family,splicing_score
colnames(vcf)<-c("chr" ,"pos","ref","alt", "ind","gene","ensg","enst","maf", "cadd_phred","impactful","impact_score","dna_change","aa_change","disease","gene_family","regional_score","num_rv" )
print(head(vcf))
vcf=vcf %>% mutate(varcat=ifelse(nchar(vcf$ref)>1 | nchar(vcf$alt) >1,"INDEL","SNP" ))
vcf_red<-vcf[,c("ind","gene","ensg","maf","cadd_phred","regional_score","varcat")]
spl<-fread(spl_f,fill=TRUE)
colnames(spl)<-c("chr","start","end","cluster","sample","z","rank","ensg","enst","ensg_withperiod", "hgnc"  ,  "genic_score" ,"V13")
spl_red<-spl[,c("sample","z","ensg","enst","genic_score")]
spl_vcf<-merge(spl_red,vcf_red,all=T,by="ensg")
spl_vcf$status<-status
spl_vcf$tissue<-tiss
print(head(spl_vcf))
write.table(spl_vcf,file=out,col.names =F,sep="\t",quote=F,row.names=F)
