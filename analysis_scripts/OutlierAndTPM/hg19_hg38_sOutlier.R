print(0)
library(data.table)
library(stringi) #chosen for supposed speediness
library(ggplot2)
library(tidyverse)
library(optparse)

print(.1)

option_list = list(
  make_option(c( "--hg19"), type="character", default=NULL, help="z hg19", metavar="character"),
  make_option(c("--hg38"), type="character", default=NULL, help="z hg38", metavar="character"),
  make_option(c("--use_nonensembl"), type="logical", default=FALSE, help="z hg38", metavar="character"),
  make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
z_file_hg19<-as.character(opt$hg19)
z_file_hg38<-as.character(opt$hg38)
use_nonsensembl<-as.logical(opt$use_nonensembl)
outfile<-as.character(opt$outfile)
print(.2)
###z values
# # 
# z_file_hg19="/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg19.sorted_dedupOptical_minMQ255.bed.gz "
# z_file_hg38="/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg38.sorted_dedupOptical_minMQ255.bed.gz"
# outfile="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/zscores_genesdiff_hg19_hg38.soutliers.txt.gz"


colnames_19cols=c("chr","start","end","cluster","sample_id","z",
                  "junc_name","score","strand","splice_site","acceptors_skipped","exons_skipped",
                  "donors_skipped","anchor","known_donor","known_junction","transcript","gene","x")
print(.3)
colnames_11cols=c("chr","start","end","cluster","sample_id","z",
                  "transcript","gene","chm13_gene","collapsedon","x")
colnames_9cols=c("chr","start","end","cluster","sample_id","z",
                  "transcript","gene","x")
print(0.4)
 z_hg19<-fread(z_file_hg19)
 print(1)
 if(ncol(z_hg19)==19){
   colnames(z_hg19)<-colnames_19cols
 }else if(ncol(z_hg19)==11){
   colnames(z_hg19)<-colnames_11cols
   
 }else if(ncol(z_hg19)==9){
   colnames(z_hg19)<-colnames_9cols
   
 }else{
   stop("ERROR NOT PROPER FORMAT FOR COLS")
 }
 #colnames(z_hg19)<-c("sample_id","gene","z","rank_all","rank_underexp","rank_overexp")
 print(2)
 
 z_hg38<-fread(z_file_hg38)
 print(3)
 #colnames(z_hg38)<-c("sample_id","gene","z","rank_all","rank_underexp","rank_overexp")
 if(ncol(z_hg38)==19){
   colnames(z_hg38)<-colnames_19cols
 }else if(ncol(z_hg38)==11){
   colnames(z_hg38)<-colnames_11cols
   
 }else{
   stop("ERROR NOT PROPER FORMAT FOR COLS")
 }
  #print(head(z_hg19))
  print(head(z_hg38))
hg19<-z_hg19[,c("sample_id","gene","transcript","z","chr")] #%>%mutate("build"="hg19")
if(use_nonsensembl){
  hg38<-z_hg38[,c("sample_id","chm13_gene","transcript","z","chr")] #%>%mutate("build"="hg38")
  colnames(hg38)[2]<-"gene"
  print(colnames(hg38))

}else{
  hg38<-z_hg38[,c("sample_id","gene","transcript","z","chr")] #%>%mutate("build"="hg38")

}
print(head(hg38$gene))
print("hg 38 gene id ^^, hg 19 gene id below")
print(head(hg19$gene))
hg19$gene<-sapply(str_split(hg19$gene,"\\."),"[[",1)
hg38$gene<-sapply(str_split(hg38$gene,"\\."),"[[",1)
print(head(hg38$gene))
print("hg 38 gene id ^^, hg 19 gene id below")
print(head(hg19$gene))
print(4)
print("redone")
z_all<-merge(hg19, hg38, all=T,by=c("sample_id","gene"))
print(5)
z_all=z_all[,-c("chr.y")]
z_all$z_diff<-z_all$z.y-z_all$z.x
# z_all$rank_ind_diff<-z_all$rank_ind_hg38-z_all$rank_ind_hg19
print(head(z_all))
fwrite(z_all,outfile)
