library(data.table)
library(sva)
library(tidyverse)
# library(tximport)
library(dplyr)
# library(readr)
# library(stringr)
# library(splines)
# library(peer)
# library(plyr)
library(optparse)
# library(reshape2)
print("TPM_to_Zscore.R")
predictor_sig <- function(sva_object, predictor, cutoff) {
  #print(lm(x~factor(predictor)))
  return(which(apply(sva_object, 2, function(x)
    any(summary(lm(x ~ factor(predictor)))[[4]][-1,4] < cutoff))))
  # return(which(sapply(sva_object, function(x) 
  #   any(summary(lm(x ~ factor(predictor)))[[4]][-1,4] < cutoff)))) 
}
#--- OPTION PARSER
option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="metadata file -- UDN", metavar="character"),
  make_option(c("-t", "--tissue"), type="character", default=NULL, help="tissue type (ex: blood, fibroblast, etc) ", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-b", "--build"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-s", "--sva_fit_file"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-r", "--sva_resids_file"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-l", "--feature_level"), type="character", default=NULL, help="enst OR ensg ONLY", metavar="character"),
  make_option(c("-g", "--outfile_transformed"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-u", "--counts"), type="character", default=NULL, help="counts udn", metavar="character")
  
) 
#I personally just like saving the parameters pass through for ease 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
metadata_file<-as.character(opt$metadata)
tissue_type<-(as.character(opt$tissue))
outfile<-as.character(opt$outfile)
build<-as.character(opt$build)
outfile_transformed<-as.character(opt$outfile_transformed)
counts_file<-as.character(opt$counts)
sva_fit_file<-as.character(opt$sva_fit_file)
sva_resids_file<-as.character(opt$sva_resids_file)
#
# I just have this here so I can load things into my local Rstudio and play around
#  counts_file<-"/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results"
# tissue_type<-"Blood"
# metadata_file<-"/Volumes/groups/smontgom/shared/UDN/Data/SampleTables/sample_table_UDN_aug2022.txt"

## filter genes read count filter and tpm filter
#default values: people_filt and tpm_filt are used such that 20% of individuals must have a TPM value>0.15
#default value: min_sample_size means that if 20% of individuals is less than 30, instead 30 individuals must have TPM values > 0.5
#returns TPMs that passes these filters
filter_tpms = function(tpm, people_filt = .2, tpm_filt = 0.15,min_sample_size=30) {
  genes = tpm$ensg
  tpm_both = tpm[,-1]  %>%  mutate_all(as.numeric) #change to numeric
  ind_filt_both = round(people_filt*ncol(tpm_both)) #this is 20% of people
  if(ind_filt_both<min_sample_size){ind_filt_both=min_sample_size} #makes sure that the number of samples with a certain TPM value is at least the min_sample_size
  #how many people pass 20% of tpm filtering and min reads numbering
  #or at least 30 people
  indices_keep_both = (rowSums(tpm_both > tpm_filt) >= ind_filt_both ) #which rows pass this filter
  tpm_out_both = data.frame(tpm_both[indices_keep_both, ]) #subset TPM values to just those that pass filters
  rownames(tpm_out_both) = genes[indices_keep_both] #make sure to include gene names
  return(tpm_out_both)
}


######## Load sample info. 
# get metadata that is just for this tissue
metadata <- fread(metadata_file) %>% filter(TISSUE==tissue_type) #replace by tissue of interest
print("metadata")
#read in counts file and reformat
counts_melt<-fread(counts_file)
colnames(counts_melt)<-c("sample","ensg","enst","TPM")

##remove samples that we cant do batch correction on
#must be greater than or equal to min_samples_per_batch
min_samples_per_batch=3
samples<-unique(counts_melt$sample)
metadata_filt=metadata%>%dplyr::filter(SAMPLE %in% samples)
not_enough_samples_in_batch=names(which(table(metadata_filt$RUN)<min_samples_per_batch))
samples_to_remove=metadata_filt%>%dplyr::filter(RUN%in%not_enough_samples_in_batch)%>%select(SAMPLE) %>%unlist()

print("counts")
print("0")
counts_melt_filt<-counts_melt%>%dplyr::filter(ensg!="None"& !(sample %in%samples_to_remove))
print("1")
counts_melt_filt_uniq=counts_melt_filt #%>%group_by(sample,ensg)%>%mutate(ensg=make.unique(ensg,sep="-"))
print("2")
counts<-data.table::dcast(counts_melt_filt_uniq[,c("ensg","sample","TPM")],ensg~sample)
print("before filt")
print((counts)[1:4,1:4])
#### filter  counts to those that pass tpm and ind filter
all_counts_filtered<-filter_tpms(counts) #custom function, returns TPMs that passes these filters
all_counts_filtered_rns<-cbind(rownames(all_counts_filtered),all_counts_filtered)
colnames(all_counts_filtered_rns)[1]<-"ensg"
print("counts filtered")
print(head(all_counts_filtered_rns[1:4,1:4]))
# Moderated log transform --> so laure log10(x+1) and nicole does log2(x + 2)
dat_filter_log<-t(apply((all_counts_filtered_rns)[,-1], 1, function(x) log10(x+1)))
print("transformed")
print(head(dat_filter_log))
## Unit variance and center for cols with var != 0
dat_filter_log_rmnovar <- dat_filter_log[, apply(dat_filter_log, 2, function(x) !var(x) == 0)] 
dat_filter_logscaled_rmnovar<-t(apply(dat_filter_log_rmnovar,1,function(x){scale(x,center=T,scale=T)}))
colnames(dat_filter_logscaled_rmnovar)<-colnames(dat_filter_log_rmnovar)
# for (i in 1:ncol(dat_filter_log)) {
#   if (var(dat_filter_log[ ,i]) != 0) dat_filter_log[, i] <- scale(dat_filter_log[, i], center=TRUE, scale=TRUE)
# }
# ###scale(myresids_num, center=TRUE, scale=TRUE)
# remove genes with zero variance

dat_filter_logscaled_rmnovar_t <- t(dat_filter_logscaled_rmnovar) #  transpose
#save this this, as this is the non-corrected count data that is transformed
# write.table(dat_filter_logscaled_rmnovar_t,paste0("/Volumes/groups/smontgom/shared/UDN/PreprocessingHG19/eOutliers/transformedcounts_beforeSVA_",tissue_type,"_hg19.txt"),
#             quote = F,sep="\t", row.names=T, col.names=T)
write.table(dat_filter_logscaled_rmnovar_t, outfile_transformed, sep="\t", row.names=T, col.names=T, quote=F)
## SVA
print("putting it into this model...")
print(head(dat_filter_logscaled_rmnovar))
#prep matrix for SVA

mod <- model.matrix(~1, data=as.data.frame((dat_filter_logscaled_rmnovar_t)))
print("about to run SVA.....")
print(head(mod))
#this is the line that actually runs SVA!
sva_fit <- sva(dat_filter_log_rmnovar, mod,mod0=NULL, method="two-step") 
print(head(sva_fit))
#add back  in row/colnames to SVA
colnames(sva_fit$sv)<-paste0("SV",c(1:ncol(sva_fit$sv)))
rownames(sva_fit$sv)<-colnames(dat_filter_logscaled_rmnovar)
print("SVA FIT COMPLETED")
#save the SVA data for analyses plots
#saveRDS(sva_fit,sva_fit_file)
#paste0("/Volumes/groups/smontgom/shared/UDN/PreprocessingHG19/eOutliers/sva_fit_",tissue_type,"_",build,"RData")
sva_fit=readRDS(paste0("/Volumes/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Blood.hg19.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.svafit.RData"))

print("Regress out svs")
## Regress out SVs

metadata_filt=metadata%>%dplyr::filter(SAMPLE %in% rownames(sva_fit$sv))
metadata_filt_reorder<-metadata_filt[match(rownames(sva_fit$sv),metadata_filt$SAMPLE)]
metadata_filt_reorder=metadata_filt_reorder%>% group_by(RUN)%>%mutate(RIN = ifelse(is.na(RIN), median(RIN, na.rm=T), RIN))
metadata_filt_reorder$RUN<-paste0("RUN",metadata_filt_reorder$RUN)
with_batch=model.matrix(~RUN+RIN+sex,metadata_filt_reorder)
rownames(with_batch)=rownames(mod)
modsv <- cbind(with_batch,sva_fit$sv)

#  modsv = model.matrix(~1, data=as.data.frame(dat_filter_logscaled_rmnovar_t))
# i think this was written this way so that if you want to protect variables in SVA you can
print('fit')
# print("fit completed")
# paste0("/Volumes/groups/smontgom/shared/UDN/PreprocessingHG19/eOutliers/sva_resids_",tissue_type,"_hg19.txt"
myzs=data.frame()
for(i in 1:ncol(sva_fit$sv)){
  mycols=1:ncol(sva_fit$sv)
  mycols=mycols[mycols!=i]
  modsv <- as.matrix(sva_fit$sv[,mycols])
  fitsv <- lm.fit(modsv, dat_filter_logscaled_rmnovar_t)
  
  ## Drop rownames and make matrix numeric
  myresids=cbind(rownames(fitsv$residuals), fitsv$residuals)
  myresids_mat <- as.matrix(myresids[, -1]) # convert dataframe to matrix
  myresids_num <- apply(myresids_mat, 2, as.numeric) # convert to numeric
  rownames(myresids_num) <-unlist(myresids[, 1]) # add back rownames
  myresids_num[sample,str_detect(colnames(myresids_num),gene)]
  
  print("scaling")
  ## Scale and center
  #Part where z score is actually calculated
  gene_zscore_scale <- scale(myresids_num, center=TRUE, scale=TRUE)
  print(paste0("SV",i,": ",zscore))
  myzs<-rbind(myzs,data.frame(paste0("SV",i),zscore))
  
}
  
 
write.table( fitsv$residuals,
             sva_resids_file,
             quote = F,sep="\t", row.names=T, col.names=T)

#think about each batch
sample_to_run=metadata_filt_reorder$RUN
names(sample_to_run)=metadata_filt_reorder$SAMPLE
run_to_sample=metadata_filt_reorder$SAMPLE
names(run_to_sample)=metadata_filt_reorder$RUN
j=0
for(gene_name in colnames(gene_zscore_scale)){
  this_gene=as.data.frame(gene_zscore_scale[,gene_name])
  #print(head(this_gene))
  split_by_batch=split(this_gene,sample_to_run[rownames(this_gene)])
  batch_name=0
  for(df in split_by_batch){
    batch_name=batch_name+1
    this_batch_name=names(split_by_batch)[batch_name]
    batch_mean=mean(df[,1])
    batch_var=var(df[,1])
    myresids_num[,gene_name]
    if(batch_var<0.1){  
      #print(batch_var)
      j=j+1
      #print("HIIII");
      # print(gene_zscore_scale[,gene_name])
      samples_in_this_batch=rownames(df) #which are the samples we are excluding
      fixing_gene=myresids_num[,gene_name] #get the column we are fixing
      turn_na=rep(NA,length(samples_in_this_batch)) #create dic for NA samples to add at the fixed_gene= step
      names(turn_na)<-samples_in_this_batch
      pre_z=t(as.data.frame(fixing_gene))[,!(rownames(myresids_num)%in%samples_in_this_batch)] #get samples that don't have bad variance/mean
      z_subsetted=scale(pre_z, center=TRUE, scale=TRUE) #rescale excluding samples with bad variance/mean
      fixed_gene=as.data.frame(c(z_subsetted[,1],turn_na)) #combine fix and na samples
      fixed_gene_reorder=fixed_gene[order(match(rownames(fixed_gene), rownames(myresids_num))), , drop = FALSE] #make sure order is correct
      gene_zscore_scale[,gene_name]=fixed_gene_reorder[,1] #replace in place
    }
    #print(paste0(this_batch_name, ": ", " mean=",batch_mean, ", var=",batch_var))
  }
  #break
}


## Melt zscore matrix (reformat)
gene_zscore_scale.m=reshape2::melt(gene_zscore_scale)
colnames(gene_zscore_scale.m)=c("sample_id", "gene","zscore")
#calculate under/over expression ranks (genetic counselors requested this!)
zscores_ranked<-gene_zscore_scale.m %>% mutate(rank_all=rank(-abs(zscore))) %>% mutate(rank_underexp=rank(zscore)) %>% mutate(rank_overexp=rank(-zscore))

## Write data
print(paste0("writing to: ",outfile))
write.table(zscores_ranked, outfile, sep="\t", row.names=F, col.names=T, quote=F) 
