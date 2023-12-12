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
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-g", "--outfile_transformed"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-c", "--counts_control"), type="character", default=NULL, help="counts control", metavar="character"),
  make_option(c("-u", "--counts_udn"), type="character", default=NULL, help="counts udn", metavar="character")

) 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
metadata_file_udn<-as.character(opt$metadata_udn)
metadata_file_control<-as.character(opt$metadata_control)
tissue_type<-(as.character(opt$tissue))
# output_file_spline<-as.character(opt$out_spline)
outfile<-as.character(opt$outfile)
outfile_transformed<-as.character(opt$outfile_transformed)
counts_udn_file<-as.character(opt$counts_udn)
counts_control_file<-as.character(opt$counts_control)

# 
# metadata_file_udn<-"/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/2019_12_05_Rare_Disease_Metadata.tsv"
 # metadata_file_control<-"/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
# metadata_file_control_phenos<-"/Volumes/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
# tissue_type<-"Fibroblast"
# counts_udn_file<-"/Volumes/groups/smontgom/shared/UDN/PreprocessingHG38/eOutliers/Fibroblast_combined.genes.results"
# counts_control_file<-"/Volumes/groups/smontgom/shared/UDN/PreprocessingHG38/eOutliers/gtex_Fibroblast_hg38.tpm.txt.gz"

ztrans.filter = function(tpm, people.filt = .2, tpm.filt = 0.15) {
  # genes_vec = tpm$ensg
  #Both
  tpm.both = tpm[,-1]  %>%  mutate_all(as.numeric) #[, c(covs.both), with = F]
  ind.filt.both = round(people.filt*ncol(tpm.both)) #20% of people
  if(ind.filt.both<30){ind.filt.both=30}
  #how many people pass 20% of tpm filtering and min reads numbering
  #or at least 30 people
  #this passes female or male
  indices.keep.both = (rowSums(tpm.both > tpm.filt) >= ind.filt.both ) 
  tpm.out.both = data.frame(tpm.both[indices.keep.both, ])
  ####tpm.out.both = scale(t(log2(tpm.cut.both + 2))) #log and z transform
  rownames(tpm.out.both) = as.character(tpm$ensg[which(indices.keep.both==TRUE)])
  return(tpm.out.both)
}


######## Load sample info. 
#### Filters metadata
print("load udn")
metadata_udn <- read_tsv(metadata_file_udn) %>% filter(institution!="CHEO",status!="FAILED") %>% filter(source_of_RNA==tissue_type) #replace by tissue of interest
head(metadata_udn)
sex_dic=metadata_udn$sex
names(sex_dic)=metadata_udn$sample_id
#get only whole blood and rna batch details for metadata
# metadata_control <- read_tsv(metadata_file_control) %>% filter(SMTSD == "Whole Blood" )  %>% filter(grepl("RNA",SMNABTCHT)) %>% as.data.frame() #gtex batch: SMNABTCH, center=SMCENTER
#### Filters counts and reformat
counts_udn_melt<-fread(counts_udn_file) #######FILTER SO HAS SAME NAMES AS METADATA
samples_in_counts<-unique(counts_udn_melt$sample)
samples_in_counts[grepl("SC",colnames(counts_udn_melt))]<-paste0("RD0",as.numeric(str_sub(samples_in_counts[grepl("SC",colnames(counts_udn_melt))],-2,-1))+12)
failed_ids=samples_in_counts[which(samples_in_counts %in% metadata_udn$sample_id==FALSE)]
counts_udn_melt_rmFAILS<-counts_udn_melt%>% dplyr::filter(!sample %in% failed_ids)
counts_udn<-data.table::dcast(counts_udn_melt_rmFAILS[,c("ensg","sample","TPM")],ensg~sample)
metadata_udn<-metadata_udn %>% dplyr::filter(sample_id %in% sort(samples_in_counts) )
sample_id_udn <-colnames(counts_udn)[-1] #########subset to just name and TPM
# counts_udn$ensg<-sapply(strsplit((counts_udn$ensg),"\\."),"[[",1)
#deal w things like ensg_par_y
split_ensgs<-strsplit(counts_udn%>% pull(ensg),"_")
counts_udn$ensg=unlist(lapply(split_ensgs,function(y){if(length(y)==3){paste0(sapply(strsplit(unlist(y)[1],"\\."),"[[",1),"_",unlist(y)[2],"_",unlist(y)[3])}else{sapply(strsplit(y,"\\."),"[[",1)}}))
#colnames(counts_udn)<-c(samples_in_counts)
# counts_udn<-counts_udn[,c("ensg",metadata_udn$sample_id)]
print("load control")

#get only whole blood and rna batch details for metadata
tissue_type_dic<-c("Whole Blood","Muscle - Skeletal","Cells - Cultured fibroblasts","Heart - Left Ventricle")
 names(tissue_type_dic)<- c("Blood","Muscle","Fibroblast","Heart")
metadata_control <- fread(metadata_file_control) %>% filter(SMTSD == tissue_type_dic[tissue_type])  %>% filter(grepl("RNA",SMNABTCHT)) %>% as.data.frame() #gtex batch: SMNABTCH, center=SMCENTER
metadata_pheno_control<-fread(metadata_file_control_phenos)
sexdic_control=ifelse(metadata_pheno_control$SEX==2,"F","M")
names(sexdic_control)<-gsub("-",".",metadata_pheno_control$SUBJID)
sexdic_all<-c(sexdic_control,sex_dic)
counts_control<-fread(counts_control_file)
colnames(counts_control)[1]<-"ensg"
# counts_control$ensg<-sapply(strsplit((counts_control$ensg),"\\."),"[[",1)
split_ensgs_control<-strsplit(counts_control%>% pull(ensg),"_")
counts_control$ensg=unlist(lapply(split_ensgs_control,function(y){if(length(y)==3){paste0(sapply(strsplit(unlist(y)[1],"\\."),"[[",1),"_",unlist(y)[2],"_",unlist(y)[3])}else{sapply(strsplit(y,"\\."),"[[",1)}}))

#counts_control$ensg<-sapply(strsplit((counts_control$ensg),"\\."),"[[",1)
sample_id_control<-colnames(counts_control)[-1]
#counts_control_melt<-melt(counts_control) #change into format similar to other one
# colnames(counts_control)<-c("ensg","rd_id","TPM")
#### correct naming, subset metadata to just the ids necessary for the future <3 
metadata_control$SAMPID<-sapply(strsplit(metadata_control$SAMPID,"-"),function(x){paste0(x[1],"-",x[2])})
metadata_control<-metadata_control %>% dplyr::filter(SAMPID %in% colnames(counts_control))
print("COMBINING")
print(head(counts_control))
print("UDN")
print(head(counts_udn))
#### merge the counts together
all_counts<-merge(counts_udn,counts_control,by=c("ensg"),all=T)
sample_id<-c(sample_id_udn,sample_id_control)
print(head(all_counts))
#all_counts<-rbind(cbind(counts_udn,ind_status="case"),cbind(counts_control,ind_status="control))
# all_counts<-rbind(counts_udn[,c("ensg","rd_id","TPM")],counts_control) 
# colnames(all_counts)<-c("ensg","sample_id","TPM")
# all_counts<-rbind(counts_udn,counts_control)

#### filter  counts to those that pass tpm and ind filter
all_counts_filtered<-ztrans.filter(all_counts)
all_counts_filtered_rns<-cbind(rownames(all_counts_filtered),all_counts_filtered)
colnames(all_counts_filtered_rns)[1]<-"ensg"

# Moderated log transform --> so laure log10(x+1) and nicole does log2(x + 2)
#dat_filter<-melt(all_counts_filtered_rns,id="ensg")
dat_filter_log<-apply((all_counts_filtered_rns)[,-1], 1, function(x) log10(x+1))
#dat_filter_log <- dat_filter[, apply(all_counts_filtered_rns[,-1], 1, function(x) log10(x+1))]
# dat_filter_log <- cbind(dat_filter[, 1], dat_filter_log)

## Unit variance and center for cols with var != 0
temp_dat <- t(dat_filter_log[ ,-1])
#rownames(temp_dat) <- as.character((all_counts_filtered_rns[,1]))

#transform to normal if there is variance
for (i in 1:ncol(temp_dat)) {
	if (var(temp_dat[ ,i]) != 0) temp_dat[, i] <- scale(temp_dat[, i], center=TRUE, scale=TRUE)
}

temp_dat_var <- temp_dat[, apply(temp_dat, 2, function(x) !var(x) == 0)] # remove genes with zero variance

dat_filter_log_scale <- (t(temp_dat_var)) # rename and transpose
write.table(dat_filter_log_scale, outfile_transformed, sep="\t", row.names=T, col.names=T, quote=F) 

print("SVA")

## SVA
print("putting it into this model...")
mod <- model.matrix(~1, data=as.data.frame((dat_filter_log_scale)))
print("about to run SVA.....")
sva_fit <- sva(t(dat_filter_log_scale), mod, method="two-step") ##TRANSPOSE?
print("SVA FIT COMPLETED")
saveRDS(sva_fit,"/oak/stanford/groups/smontgom/shared/UDN/Preprocessing/eOutliers/sva_fit_WholeBlood.RData")
## Splines
# ################## ADDD IN GTEX BATCH OKAY!!! GREAT :) 
# batch <- c(metadata_udn$batch, metadata_control$SMNABTCH) #gtex batch: SMNABTCH, center=SMCENTER
# study <- c(metadata_udn$institution, metadata_control$SMCENTER)
# 
# batch_pred_sig<-predictor_sig(sva_fit$sv, droplevels(factor(batch)), 1e-30)
# study_pred_sig<-predictor_sig(sva_fit$sv, droplevels(factor(study)), 1e-30)
# sig_sv <- unique(c(batch_pred_sig,study_pred_sig))


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

# fitsv<-do_fit(sig_sv,sva_fit,modsv,dat_filter_log_scale,spline=TRUE)
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

# gene_zscore_scale_spline<-get_zscore(fitsv)
gene_zscore_scale_nospline<-get_zscore(fitsv_nospline)

## Write data
print(paste0("writing to: ",outfile," and ", outfile))
#write.table(gene_zscore_scale.m, "2019_12_05_Muscle_outliers_zscore_pair_spline.txt", sep="\t", row.names=F, col.names=T, quote=F) 
write.table(gene_zscore_scale_nospline, outfile, sep="\t", row.names=F, col.names=T, quote=F) 
# write.table(gene_zscore_scale_nospline, output_file_nospline, sep="\t", row.names=F, col.names=T, quote=F) 


