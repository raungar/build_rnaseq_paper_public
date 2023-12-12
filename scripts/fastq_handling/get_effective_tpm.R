library(data.table)
library(dplyr)
library(optparse)

#--- OPTION PARSER
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, help="infile", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="outfile", metavar="character")
) 
#I personally just like saving the parameters pass through for ease 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


rsem<-fread(opt$infile)%>%dplyr::filter(effective_length>0)
across_sample_scalingfactor<-rsem%>%mutate(effective_length_kb = effective_length/1000) %>%
  mutate(expected_rpk = expected_count/effective_length_kb)%>%
  pull(expected_rpk)%>%sum()/10^6
rsem_expectedtpm = rsem %>%
  mutate(expected_rpk = expected_count/effective_length/1000) %>%
  mutate(expected_tpm = round(expected_rpk/across_sample_scalingfactor,2)) %>%
  select(-expected_rpk)

fwrite(rsem_expectedtpm,file=opt$outfile,sep="\t",quote=F)


