library(data.table)
library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("--infile"), type="character", default=NULL, help="infile", metavar="character"),
  make_option(c("--pull_type"), type="character", default=NULL, help="enst or ensg", metavar="character"),
  make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile=opt$infile
pull_type=opt$pull_type
outfile=opt$outfile
print(pull_type)

data=fread(infile)
colnames(data)<-c("chr","start","end","clus","RD","z","rank","ensg","enst")
head(data)
collapsed=data%>% group_by_at(pull_type)  %>% slice(which.max(z))
write_csv(collapsed,file=outfile,quote=F)
