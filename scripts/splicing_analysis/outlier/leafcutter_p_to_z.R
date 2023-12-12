library(data.table)
library(GenomicRanges)
library(optparse)
library("tidyverse")



option_list = list(
  make_option(c( "--infile"), type="character", default=NULL, help="udn_blood_pVals.txt", metavar="character"),
  make_option(c( "--outfile"), type="character", default=NULL, help="udn_blood_pVals.txt", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
print("next")
opt = parse_args(opt_parser)
infile=as.character(opt$infile)
outfile=as.character(opt$outfile)


#read in actual infile
#infile<-"/Volumes/groups/smontgom/shared/UDN/PreprocessingHG19/Splicing/LeafCutterMD/udn_Fibroblast_pVals.txt"
this_header<-fread(cmd=paste('head -1 ', infile))
pvals<-fread(infile,header=F)
colnames(pvals)<-c("pos",colnames(this_header))
print(pvals[1:10,1:10])

#### padj
#adjust the pvals
number_of_cluster<-length(unique(sapply(strsplit(pvals$pos,":"),"[[",4)))
number_of_samples<-ncol(pvals)-1
number_of_tests<-number_of_samples*number_of_cluster
adjusted<-lapply((pvals[,2:ncol(pvals)]),function(x){p.adjust(x,method="fdr",n=number_of_tests)})
pvals_adjusted<-cbind(pvals[,1],do.call(cbind.data.frame,(adjusted)))
print(pvals_adjusted[1:5,1:14])

##z score
zscores<-cbind(pvals[,1],
                do.call(cbind.data.frame,
                        lapply((pvals_adjusted[,2:ncol(pvals_adjusted)]),function(x){abs(qnorm(x/2))})
                        )
             )
zscore_melted<-melt(zscores)
colnames(zscore_melted)<-c("pos","sample","zscore")
print("z score")
print(zscore_melted[1:10,])
long_data<-zscore_melted 
print(as.data.frame(head(long_data,40)))
fwrite(long_data,outfile)