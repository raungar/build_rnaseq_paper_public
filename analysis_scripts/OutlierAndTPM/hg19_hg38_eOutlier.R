library(data.table)
library(stringi) #chosen for supposed speediness
library(ggplot2)
library(tidyverse)
library(optparse)



option_list = list(
  make_option(c( "--build1"), type="character", default=NULL, help="z build1", metavar="character"),
  make_option(c("--build2"), type="character", default=NULL, help="z build2", metavar="character"),
    make_option(c( "--build1_name"), type="character", default=NULL, help="z build1", metavar="character"),
  make_option(c("--build2_name"), type="character", default=NULL, help="z build2", metavar="character"),
  make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
z_file_build1<-as.character(opt$build1)
z_file_build2<-as.character(opt$build2)
build1_name<-as.character(opt$build1_name)
build2_name<-as.character(opt$build2_name)
outfile<-as.character(opt$outfile)

###z values
# 
# z_file_build1="/Volumes/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Blood.hg19.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt"
# z_file_build2="/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt"
#build1_name="hg19";build2_name="hg38";
#outfile="/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz"

 z_build1<-fread(z_file_build1)
 colnames(z_build1)<-c("sample_id","gene","z","rank_all","rank_underexp","rank_overexp")
 z_build2<-fread(z_file_build2)
 colnames(z_build2)<-c("sample_id","gene","z","rank_all","rank_underexp","rank_overexp")
  head(z_build1)
  head(z_build2)
build1<-z_build1[,c("sample_id","gene","z")] #%>%mutate("build"="build1")
build2<-z_build2[,c("sample_id","gene","z")] #%>%mutate("build"="build2")
build2_repeats=build2%>%dplyr::filter(str_detect(gene,"-"))%>%mutate(gene_orig=str_replace(gene,"-[0-9]*",""))

build2=build2%>%dplyr::filter(!str_detect(gene,"-"))
build2$gene<-sapply(str_split((build2$gene),"_"),"[[",1) ##str_replace(build2$gene,"_[0-9]+","")
build1$gene<-str_replace(build1$gene,"_[0-9]+","") # sapply(str_split((build1$gene),"_"),"[[",1)
z_all<-merge(build1, build2, all=T,by=c("sample_id","gene"))
z_repeats<-merge(build1,build2_repeats,by.x=c("sample_id","gene"),by.y=c("sample_id","gene_orig"),all.y=T)
z_repeats_red<-z_repeats[,c(1,4,3,5)]
colnames(z_repeats_red)<-colnames(z_all)
z_all<-rbind(z_all,z_repeats_red)
z_all<-z_all%>%mutate(z_diff=z.y - z.x)
colnames(z_all)<-c("sample_id","gene",paste0("z_",build1_name),paste0("z_",build2_name),"z_diff")

# z_all$z_diff<-z_all$z_build2-z_all$z_build1
fwrite(z_all,outfile)
