library(data.table)
library(dplyr)
library(optparse)
option_list = list(
  make_option(c( "--infile"), type="character", default=NULL, help="z hg19", metavar="character"),
  make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile<-as.character(opt$infile)
outfile<-as.character(opt$outfile)
print(date())
#infile="/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/Mapping/mapping_hg19_ID.txt"
mapping=fread(infile,header=FALSE)
print(date())
setnames(mapping,c("build","sample","enst","NH","num_mapped"))
# colnames(mapping)=c("build","sample","enst","NH","num_mapped")
print(date())
setkey(mapping, num_mapped)
print(date())

# mapping_nozero <- mapping[num_mapped > 1]
# mapping_nonzero=mapping[.( c(1,600000))]
# mapping_nozero=subset(mapping,num_mapped>0)
# mapping_nozero=mapping[which(mapping$num_mapped>0),]
print(date())
mapping[,mapping_status:=ifelse(mapping$num_mapped==1,"uniquely_mapped","multi_mapped")]
                         #ifelse(mapping$num_mapped==0,"unmapped","uniquely_mapped"))]
print(date())
mapping_summarized=mapping[, .(count = .N, uniquely_mapped = sum(mapping_status=="uniquely_mapped"),multi_mapped=sum(mapping_status=="multi_mapped")), by = c("build","sample","enst")]
print(date())
# mapping_summarized=mapping_nozero%>%
#     group_by(build,sample,enst)%>%
#     mutate(mapping_status=ifelse(num_mapped==1,"uniquely_mapped","multi_mapped"))%>%
#     summarise(total_reads=n(),uniquely_mapped=sum(mapping_status=="uniquely_mapped"),multi_mapped=sum(mapping_status=="multi_mapped"))

fwrite(mapping_summarized,file=outfile,quote=FALSE,sep="\t")
print(date())
