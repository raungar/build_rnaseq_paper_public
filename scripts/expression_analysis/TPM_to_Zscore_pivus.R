library(data.table)
library(sva)
library(tximport)
library(dplyr)
library(readr)
library(stringr)
library(splines)
library(plyr)
library(optparse)
library(reshape2)

predictor_sig <- function(sva_object, predictor, cutoff) {
	#print(lm(x~factor(predictor)))
	return(which(apply(sva_object, 2, function(x)
		any(summary(lm(x ~ factor(predictor)))[[4]][-1,4] < cutoff))))
  # return(which(sapply(sva_object, function(x) 
  #   any(summary(lm(x ~ factor(predictor)))[[4]][-1,4] < cutoff)))) 
}
#--- OPTION PARSER
option_list = list(
  make_option(c("-m", "--metadata_udn"), type="character", default=NULL, help="metadata file -- UDN", metavar="character"),
  make_option(c("-n", "--metadata_control"), type="character", default=NULL, help="metadata file -- control", metavar="character"),
  make_option(c("-t", "--tissue"), type="character", default=NULL, help="tissue type (ex: blood, fibroblast, etc) ", metavar="character"),
  make_option(c("-s", "--out_spline"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-o", "--out_nospline"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-c", "--counts_control"), type="character", default=NULL, help="counts control", metavar="character"),
  make_option(c("-u", "--counts_udn"), type="character", default=NULL, help="counts udn", metavar="character")

) 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
metadata_file_udn<-as.character(opt$metadata_udn)
metadata_file_control<-as.character(opt$metadata_control)
tissue_type<-(as.character(opt$tissue))
output_file_spline<-as.character(opt$out_spline)
output_file_nospline<-as.character(opt$out_nospline)
counts_udn_file<-as.character(opt$counts_udn)
counts_control_file<-as.character(opt$counts_control)


# metadata_file_udn<-"/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/2019_12_05_Rare_Disease_Metadata.tsv"
# metadata_file_control<-"/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/PIVUS_RNASequencingInfo.csv"
# tissue_type<-"Blood"
# counts_udn_file<-"/oak/stanford/groups/smontgom/shared/UDN/Preprocessing/eOutliers/all_udn.tpm.txt.gz"
# counts_control_file<-"/oak/stanford/groups/smontgom/shared/UDN/Preprocessing/eOutliers/all_pivus.tpm.txt.gz"

## filter genes read count filter and tpm filter
ztrans.filter = function(tpm, people.filt = .2, tpm.filt = 0.15) {
  genes = tpm$gene_id
  
  #Both
  tpm.both = tpm[,-1]  %>%  mutate_all(as.numeric) #[, c(covs.both), with = F]
  ind.filt.both = round(people.filt*ncol(tpm.both)) #20% of people
  #how many people pass 20% of tpm filtering and min reads numbering
  #this passes female or male
  indices.keep.both = (rowSums(tpm.both > tpm.filt) >= ind.filt.both ) 
  tpm.out.both = tpm.both[indices.keep.both, ]
  ####tpm.out.both = scale(t(log2(tpm.cut.both + 2))) #log and z transform
  rownames(tpm.out.both) = genes[indices.keep.both]
  return(tpm.out.both)
}


######## Load sample info. 
#### Filters metadata
metadata_udn <- read_tsv(metadata_file_udn) %>% filter(institution!="CHEO",status!="FAILED") %>% filter(source_of_RNA==tissue_type) #replace by tissue of interest
#get only whole blood and rna batch details for metadata
metadata_control <- read_csv(metadata_file_control) #batch names is in column "RUN"
#### Filters counts and reformat
counts_udn_melt<-read_tsv(counts_udn_file) #######FILTER SO HAS SAME NAMES AS METADATA
counts_melt_colnames<-colnames(counts_udn_melt)
counts_udn<-data.table::dcast(counts_udn_melt[,c("gene_id","rd_id","TPM")],gene_id~rd_id)
counts_control_melt<-read_tsv(counts_control_file,col_names = F)
colnames(counts_control_melt)<- counts_melt_colnames
counts_control_melt$gene_id<-sapply(strsplit((counts_control_melt$gene_id),"\\."),"[[",1)
counts_control<-data.table::dcast(counts_control_melt[,c("gene_id","rd_id","TPM")],gene_id~rd_id)
sample_id_control<-colnames(counts_control)[-1]
#counts_control_melt<-melt(counts_control) #change into format similar to other one
# colnames(counts_control)<-c("gene_id","rd_id","TPM")
#### correct naming, subset metadata to just the ids necessary for the future <3 
####metadata_control$SAMPID<-sapply(strsplit(metadata_control$SAMPID,"-"),function(x){paste0(x[1],"-",x[2])})
metadata_control<-metadata_control %>% dplyr::filter(RNAseq_ID %in% colnames(counts_control))
samples_in_counts<-colnames(counts_udn)
samples_in_counts[grepl("SC",colnames(counts_udn))]<-paste0("RD0",as.numeric(str_sub(samples_in_counts[grepl("SC",colnames(counts_udn))],-2,-1))+12)
metadata_udn<-metadata_udn %>% dplyr::filter(sample_id %in% sort(samples_in_counts) )
sample_id_udn <-metadata_udn$sample_id #########subset to just name and TPM
colnames(counts_udn)<-c(samples_in_counts)
counts_udn<-counts_udn[,c("gene_id",metadata_udn$sample_id)]
#### merge the counts together
all_counts<-merge(counts_udn,counts_control,by=c("gene_id"))
sample_id<-c(sample_id_udn,sample_id_control)

#all_counts<-rbind(cbind(counts_udn,ind_status="case"),cbind(counts_control,ind_status="control))
# all_counts<-rbind(counts_udn[,c("gene_id","rd_id","TPM")],counts_control) 
# colnames(all_counts)<-c("gene_id","sample_id","TPM")
# all_counts<-rbind(counts_udn,counts_control)

#### filter  counts to those that pass tpm and ind filter
all_counts_filtered<-ztrans.filter(all_counts)
all_counts_filtered_rns<-cbind(rownames(all_counts_filtered),all_counts_filtered)
colnames(all_counts_filtered_rns)[1]<-"gene_id"

# Moderated log transform --> so laure log10(x+1) and nicole does log2(x + 2)
#dat_filter<-melt(all_counts_filtered_rns,id="gene_id")
dat_filter_log<-apply((all_counts_filtered_rns)[,-1], 1, function(x) log10(x+1))
#dat_filter_log <- dat_filter[, apply(all_counts_filtered_rns[,-1], 1, function(x) log10(x+1))]
# dat_filter_log <- cbind(dat_filter[, 1], dat_filter_log)

## Unit variance and center for cols with var != 0
temp_dat <- t(dat_filter_log[ ,-1])
#rownames(temp_dat) <- as.character((all_counts_filtered_rns[,1]))

for (i in 1:ncol(temp_dat)) {
	if (var(temp_dat[ ,i]) != 0) temp_dat[, i] <- scale(temp_dat[, i], center=TRUE, scale=TRUE)
}

temp_dat_var <- temp_dat[, apply(temp_dat, 2, function(x) !var(x) == 0)] # remove genes with zero variance

dat_filter_log_scale <- (t(temp_dat_var)) # rename and transpose


## SVA
print("putting it into this model...")
mod <- model.matrix(~1, data=as.data.frame((dat_filter_log_scale)))
print("about to run SVA.....")
sva_fit <- sva(t(dat_filter_log_scale), mod, method="two-step") ##TRANSPOSE?
print("SVA FIT COMPLETED")
saveRDS(sva_fit,"/oak/stanford/groups/smontgom/shared/UDN/Preprocessing/eOutliers/sva_fit_WholeBlood_pivus.RData")
## Splines
################## ADDD IN GTEX BATCH OKAY!!! GREAT :) 
batch <- c(metadata_udn$batch, metadata_control$RUN) #gtex batch: SMNABTCH, center=SMCENTER
study <- c(metadata_udn$institution,rep("SINGLE_INSTUTION",nrow(metadata_control)))

batch_pred_sig<-predictor_sig(sva_fit$sv, droplevels(factor(batch)), 1e-30)
study_pred_sig<-predictor_sig(sva_fit$sv, droplevels(factor(study)), 1e-30)
sig_sv <- unique(c(batch_pred_sig,study_pred_sig))


print("Regress out svs")
## Regress out SVs
modsv <- cbind(mod, sva_fit$sv)
modsv_no_spline <- modsv
do_fit <- function(sig_sv,sva_fit,modsv,dat_filter_log_scale,spline){
  if(spline==TRUE){for (i in sig_sv) modsv <- cbind(modsv, bs(sva_fit$sv[, i], df=60, degree=1))}
  fitsv <- lm.fit(modsv, (dat_filter_log_scale))
  print("fit completed")
  head(fitsv)
  #save(fitsv,file="/oak/stanford/groups/smontgom/shared/UDN/Preprocessing/eOutliers/sva_fit_spline_WholeBlood.RData")
  return(fitsv)
}

fitsv<-do_fit(sig_sv,sva_fit,modsv,dat_filter_log_scale,spline=TRUE)
fitsv_nospline<-do_fit(sig_sv,sva_fit,modsv_no_spline,dat_filter_log_scale,spline=FALSE)



get_zscore<-function(fitsv){
  gene_zscore=cbind(rownames(fitsv$residuals), fitsv$residuals)
  ## Drop rownames and make matrix numeric
  gene_zscore_mat <- as.matrix(gene_zscore[, -1]) # convert dataframe to matrix
  gene_zscore_num <- apply(gene_zscore_mat, 2, as.numeric) # convert to numeric
  rownames(gene_zscore_num) <- unlist(gene_zscore[, 1]) # add back rownames
  #colnames(gene_zscore_num) <- unlist(strsplit(as.character(colnames(gene_zscore_num)), '[.]'))[c(TRUE,FALSE)] # Strip suffix from gene IDs in columns
  print("scaling")
  ## Scale and center
  gene_zscore_scale <- scale(gene_zscore_num, center=TRUE, scale=TRUE)
  
  print("melting")
  ## Melt zscore matrix
  gene_zscore_scale.m=reshape2::melt(gene_zscore_scale)
  colnames(gene_zscore_scale.m)=c("sample_id", "gene","zscore")
  return(gene_zscore_scale.m)	
}

gene_zscore_scale_spline<-get_zscore(fitsv)
gene_zscore_scale_nospline<-get_zscore(fitsv_nospline)

## Write data
print(paste0("writing to: ",output_file_spline))
#write.table(gene_zscore_scale.m, "2019_12_05_Muscle_outliers_zscore_pair_spline.txt", sep="\t", row.names=F, col.names=T, quote=F) 
#write.table(gene_zscore_scale.m, output_file, sep="\t", row.names=F, col.names=T, quote=F) 
write.table(gene_zscore_scale_spline, output_file_spline, sep="\t", row.names=F, col.names=T, quote=F) 
write.table(gene_zscore_scale_nospline, output_file_nospline, sep="\t", row.names=F, col.names=T, quote=F) 

