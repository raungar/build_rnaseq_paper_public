library(data.table)
library(stringi) #chosen for supposed speediness
library(ggplot2)
library(tidyverse)
library(optparse)



option_list = list(
  make_option(c( "--tpm_hg19"), type="character", default=NULL, help="tpm hg19", metavar="character"),
  make_option(c("--tpm_hg38"), type="character", default=NULL, help="tpm hg38", metavar="character"),
  make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
tpm_dgn_file_hg19<-as.character(opt$tpm_hg19)
tpm_dgn_file_hg38<-as.character(opt$tpm_hg38)
outfile<-as.character(opt$outfile)

 tpm_dgn_hg19<-fread(tpm_dgn_file_hg19)
 colnames(tpm_dgn_hg19)<-c("rd_id","gene_id","transcript_id","TPM")
 tpm_dgn_hg38<-fread(tpm_dgn_file_hg38)
 colnames(tpm_dgn_hg38)<-c("rd_id","gene_id","transcript_id","TPM")
  head(tpm_dgn_hg19)
  head(tpm_dgn_hg38)
dgn_hg19<-tpm_dgn_hg19[,c("rd_id","gene_id","TPM")]
dgn_hg38<-tpm_dgn_hg38[,c("rd_id","gene_id","TPM")]
head(dgn_hg38$gene_id)
print("hg 38 gene id ^^, hg 19 gene id below")
head(dgn_hg19$gene_id)
dgn_hg19$gene_id<-sapply(str_split(dgn_hg19$gene_id,"\\."),"[[",1)
dgn_hg38$gene_id<-sapply(str_split(dgn_hg38$gene_id,"\\."),"[[",1)
print("redone")
head(dgn_hg19$gene_id)
hg19_names=c("gene_id","min","iqr1","med","mean","iqr3","max","sd")
# hg38_names=c("gene_id","min_hg38","iqr1_hg38","med_hg38","mean_hg38","iqr3_hg38","max_hg38","sd_hg38")
# dgn_hg19_medgene<-dgn_hg19[,eval(hg19_names):=append(as.list(summary(TPM)[1:6]),c(sd(TPM))) ,by=c("gene_id")]
# dgn_hg38_medgene<-dgn_hg38[,eval(hg19_names):=append(as.list(summary(TPM)[1:6]),c(sd(TPM))) ,by=c("gene_id")]
dgn_hg19_medgene<-dgn_hg19[,as.list(summary(TPM)[1:6]) ,by=c("gene_id")]
dgn_hg38_medgene<-dgn_hg38[,as.list(summary(TPM)[1:6]) ,by=c("gene_id")]
print(head(dgn_hg19_medgene))
print("hg38")
print(head(dgn_hg38_medgene))

dgn_medTPM_summary<-merge(dgn_hg19_medgene, dgn_hg38_medgene, all=T,by=c("gene_id"))
print(head(dgn_medTPM_summary))

print("NOW COL")
dgn_hg19_medgene_sd<-dgn_hg19[,sd(TPM) ,by=c("gene_id")]
dgn_hg38_medgene_sd<-dgn_hg38[,sd(TPM),by=c("gene_id")]
dgn_medTPM_sd<-merge(dgn_hg19_medgene_sd, dgn_hg38_medgene_sd, all=T,by=c("gene_id"))
print(head(dgn_medTPM_sd))
dgn_medTPM_all<-cbind(dgn_medTPM_summary,dgn_medTPM_sd)


print("all")
head(dgn_medTPM_all)
colnames(dgn_medTPM_all)<-c("gene_id","min_hg19","iqr1_hg19","med_hg19","mean_hg19","iqr3_hg19","max_hg19",
					"min_hg38","iqr1_hg38","med_hg38","mean_hg38","iqr3_hg38","max_hg38","gene_id_repeat","sd_hg19", "sd_hg38")
# colnames(dgn_medTPM_all)<-c("gene_id","min_hg19","iqr1_hg19","med_hg19","mean_hg19","iqr3_hg19","max_hg19",
#                             "min_hg38","iqr1_hg38","med_hg38","mean_hg38","iqr3_hg38","max_hg38")
print("w colnames")
head(dgn_medTPM_all)
dgn_medTPM_all$iqr_range_hg19<-dgn_medTPM_all$iqr3_hg19-dgn_medTPM_all$iqr1_hg19
dgn_medTPM_all$iqr_range_hg38<-dgn_medTPM_all$iqr3_hg38-dgn_medTPM_all$iqr1_hg38
dgn_medTPM_all$iqr_range_diff<-dgn_medTPM_all$iqr_range_hg38-dgn_medTPM_all$iqr_range_hg19
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
dgn_medTPM_all$tpm_diff<-dgn_medTPM_all$med_hg38-dgn_medTPM_all$med_hg19
head(dgn_medTPM_all)
fwrite(dgn_medTPM_all,outfile)
