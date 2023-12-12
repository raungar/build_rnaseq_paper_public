library(data.table)
library(dplyr)
library(optparse)

#--- OPTION PARSER
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, help="tpm file", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="median tpm file", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile<-as.character(opt$infile)
outfile<-as.character(opt$outfile)

# infile="/Volumes/groups/smontgom/shared/UDN/PreprocessingHG38/eOutliers/gtex_Fibroblast_hg38.tpm.txt.gz"
tpm=fread(infile)

med_tpm<-apply(tpm[,-1], 1, median)
med_df=data.frame(cbind(tpm[,1],med_tpm))

fwrite(x=med_df,file=outfile,quote=F,sep="\t",row.names=F,col.names = F)