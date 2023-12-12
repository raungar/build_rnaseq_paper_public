library(data.table)
library(dplyr)
library(sva)
library(optparse)

#--- OPTION PARSER
option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="metadata file -- UDN", metavar="character"),
  make_option(c("-t", "--tissue"), type="character", default=NULL, help="tissue type (ex: blood, fibroblast, etc) ", metavar="character"),
  make_option(c("-v", "--num_svs"), type="character", default=NULL, help="tissue type (ex: blood, fibroblast, etc) ", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-b", "--build"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-r", "--resids_file"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-l", "--feature_level"), type="character", default=NULL, help="enst OR ensg ONLY", metavar="character"),
  make_option(c("-g", "--outfile_transformed"), type="character", default=NULL, help="output file", metavar="character"),
  make_option(c("-u", "--counts"), type="character", default=NULL, help="counts udn", metavar="character")
  
) 
#I personally just like saving the parameters pass through for ease 

# 
# I just have this here so I can load things into my local Rstudio and play around
#  counts_file<-"/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.dedupOptical_minMQ255_rsem.genes.results"
# tissue_type<-"Fibroblast"
# metadata_file<-"/oak/stanford/groups/smontgom/shared/UDN/Data/SampleTables/sample_table_UDN_sep2023.txt"

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

filter_counts<-function(counts_file,metadata,min_samples_per_batch=3){
  counts_melt<-fread(counts_file)
  colnames(counts_melt)<-c("sample","ensg","enst","TPM")
  
  ##remove samples that we cant do batch correction on
  #must be greater than or equal to min_samples_per_batch
  
  samples<-unique(counts_melt$sample)
  metadata_filt=metadata%>%dplyr::filter(SAMPLE %in% samples)
  not_enough_samples_in_batch=names(which(table(metadata_filt$RUN)<min_samples_per_batch))
  samples_to_remove=metadata_filt%>%dplyr::filter(RUN%in%not_enough_samples_in_batch)%>%select(SAMPLE) %>%unlist()
  
  print(samples_to_remove)
  print("counts")
  print("0")
  counts_melt_filt<-counts_melt%>%dplyr::filter(ensg!="None"& !(sample %in%samples_to_remove))
  print("1")
  counts_melt_filt_n=setDT(counts_melt_filt)[, .(n = .N), 
                          by = c("sample","ensg")]
  print(nrow(counts_melt_filt_n%>%dplyr::filter(n>1)))
  counts_melt_filt_uniq=counts_melt_filt
  #counts_melt_filt_uniq=(counts_melt_filt)%>%group_by(sample,ensg)%>%mutate(n=n())%>%dplyr::filter(n==1)%>%select(-n)
  # counts_melt_filt_uniq=counts_melt_filt%>%group_by(sample,ensg)%>%mutate(ensg=make.unique(ensg,sep="-"))
  print("2")
  counts<-data.table::dcast(counts_melt_filt_uniq[,c("ensg","sample","TPM")],ensg~sample)
  print("before filt")
  print(head(counts)[1:4,1:4])
  #### filter  counts to those that pass tpm and ind filter
  all_counts_filtered<-filter_tpms(counts) #custom function, returns TPMs that passes these filters
  all_counts_filtered_rns<-cbind(rownames(all_counts_filtered),all_counts_filtered)
  colnames(all_counts_filtered_rns)[1]<-"ensg"
  print("counts filtered")
  print(head(all_counts_filtered_rns)[1:4,1:4])
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
  return(dat_filter_logscaled_rmnovar_t)
}

filter_metadata<-function(md,dat_filter_logscaled_rmnovar_t){
  metadata_filt=md%>%dplyr::filter(SAMPLE %in% rownames(dat_filter_logscaled_rmnovar_t))
  metadata_filt_reorder<-metadata_filt[match(rownames(dat_filter_logscaled_rmnovar_t),metadata_filt$SAMPLE)]
  metadata_filt_reorder=metadata_filt_reorder%>% group_by(RUN)%>%mutate(RIN = ifelse(is.na(RIN), median(RIN, na.rm=T), RIN))
  metadata_filt_reorder$RUN<-paste0("RUN",metadata_filt_reorder$RUN)
  return(metadata_filt_reorder)
}

regress_batch<-function(metadata_filt_reorder,dat_filter_logscaled_rmnovar_t,resids_file){
  with_batch=model.matrix(~RUN+RIN+sex,metadata_filt_reorder)
  rownames(with_batch)=rownames(dat_filter_logscaled_rmnovar_t)
  modsv <- with_batch
  print("fit")
  fitsv <- lm.fit(modsv, dat_filter_logscaled_rmnovar_t)
  
  write.table( fitsv$residuals,
               resids_file,
               quote = F,sep="\t", row.names=T, col.names=T)
  return(fitsv)
}

get_z_scores<-function(fitsv){
  ## Drop rownames and make matrix numeric
  myresids=cbind(rownames(fitsv$residuals), fitsv$residuals)
  myresids_mat <- as.matrix(myresids[, -1]) # convert dataframe to matrix
  myresids_num <- apply(myresids_mat, 2, as.numeric) # convert to numeric
  rownames(myresids_num) <-unlist(myresids[, 1]) # add back rownames
  
  print("scaling")
  ## Scale and center
  #Part where z score is actually calculated
  gene_zscore_scale <- scale(myresids_num, center=TRUE, scale=TRUE)
  colnames(gene_zscore_scale)<-gsub("\\..*",'',(colnames(gene_zscore_scale)))
  return(list(myresids_num,gene_zscore_scale))
}

reformat_and_write<-function(gene_zscore_scale,outfile){
  print('h')
  ## Melt zscore matrix (reformat)
  gene_zscore_scale.m=reshape2::melt(gene_zscore_scale)
  print('i')
  colnames(gene_zscore_scale.m)=c("sample_id", "gene","zscore")
  print('j')
  #calculate under/over expression ranks (genetic counselors requested this!)
  zscores_ranked<-gene_zscore_scale.m %>% mutate(rank_all=rank(-abs(zscore))) %>% mutate(rank_underexp=rank(zscore)) %>% mutate(rank_overexp=rank(-zscore))
  
  ## Write data
  print('k')
  print(paste0("writing to: ",outfile))
  #fwrite(zscores_ranked,file="/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/TMP.TMP", sep="\t", row.names=F, col.names=T, quote=F)
  # write.table(zscores_ranked,file=outfile, sep="\t", row.names=F, col.names=T, quote=F) 
  fwrite(zscores_ranked,file=outfile, sep="\t", row.names=F, col.names=T, quote=F) 
  print('l')
}

remove_low_variance<-function(metadata_filt_reorder,gene_zscore_scale,myresids_num){
  colnames(myresids_num)=gsub("\\..*",'',(colnames(myresids_num)))
  #think about each batch
  sample_to_run=metadata_filt_reorder$RUN
  names(sample_to_run)=metadata_filt_reorder$SAMPLE
  run_to_sample=metadata_filt_reorder$SAMPLE
  names(run_to_sample)=metadata_filt_reorder$RUN
  print('b')
  j=0
  print((gene_zscore_scale[1:4,1:4]))
  print('bb')
  print(head(myresids_num[1:4,1:4]))
  #for each gene
  for(gene_name in colnames(gene_zscore_scale)){
    this_gene=as.data.frame(gene_zscore_scale[,gene_name])
    # print(head(this_gene))
    #split dataframe into $ per run
    split_by_batch=split(this_gene,sample_to_run[rownames(this_gene)])
    batch_name=0
    #for each run in the batch
    for(df in split_by_batch){
      batch_name=batch_name+1
      this_batch_name=names(split_by_batch)[batch_name]
      batch_mean=mean(df[,1])
      batch_var=var(df[,1])
      #print(batch_var)
      if(is.na(batch_var)){next}
      #myresids_num[,gene_name]
      if(batch_var<0.1){  
        print('a')
        
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
  print('c')
  return(gene_zscore_scale)
}

regress_svs<-function(metadata_filt_reorder,dat_filter_logscaled_rmnovar_t){
  print("putting it into this model...")
  #prep matrix for SVA
  mod <- model.matrix(~1, data=as.data.frame((dat_filter_logscaled_rmnovar_t)))
  print("about to run SVA.....")
  dat_filter_logscaled_rmnovar=t(dat_filter_logscaled_rmnovar_t)
  #this is the line that actually runs SVA!
  sva_fit <- sva(dat_filter_logscaled_rmnovar, mod, method="two-step") 
  print("RAN SVA")
  #add back  in row/colnames to SVA
  colnames(sva_fit$sv)<-paste0("SV",c(1:ncol(sva_fit$sv)))
  rownames(sva_fit$sv)<-colnames(dat_filter_logscaled_rmnovar)
  print("SVA FIT COMPLETED")
  #save the SVA data for analyses plots
  
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
  
  fitsv <- lm.fit(modsv, dat_filter_logscaled_rmnovar_t)
  #for (i in sig_sv) modsv <- cbind(modsv, bs(sva_fit$sv[, i], df=60, degree=1))
  # straight up the same as (dat_filter_logscaled_rmnovar_t~sva_fit$sv)
  # i think this was written this way so that if you want to protect variables in SVA you can
  fitsv <- lm.fit(modsv, (dat_filter_logscaled_rmnovar_t))
  return(list(fitsv,sva_fit))
}


regress_pcs<-function(metadata_filt_reorder,dat_filter_logscaled_rmnovar_t){
  print("putting it into this model...")
  #prep matrix for SVA
  mod <- model.matrix(~1, data=as.data.frame((dat_filter_logscaled_rmnovar_t)))
  print("about to run pca.....")
  #this is the line that actually runs SVA!
  pc_fit=prcomp(dat_filter_logscaled_rmnovar_t)
# 
#   metadata_filt_reorder=metadata_filt_reorder%>% group_by(RUN)%>%mutate(RIN = ifelse(is.na(RIN), median(RIN, na.rm=T), RIN))
#   metadata_filt_reorder$RUN<-paste0("RUN",metadata_filt_reorder$RUN)
#   with_batch=model.matrix(~RUN+RIN+sex,metadata_filt_reorder)
#   rownames(with_batch)=rownames(mod)
  # modpc <- cbind(with_batch,pc_fit$x[,1:30])
   modpc <- as.matrix(pc_fit$x[,1:30])
  
  fitpc <- lm.fit(modpc, (dat_filter_logscaled_rmnovar_t))
  return(list(fitpc,pc_fit))
}


######## Load sample info. 
# get metadata that is just for this tissue
main<-function(metadata_file,tissue_type,counts_file,outfile,outfile_transformed,build,num_svs,resids_file){
  metadata <- fread(metadata_file) %>% filter(TISSUE==tissue_type) #replace by tissue of interest
  dat_filter_logscaled_rmnovar_t<-filter_counts(counts_file,metadata,2)
  #write.table(dat_filter_logscaled_rmnovar_t, outfile_transformed, sep="\t", row.names=T, col.names=T, quote=F)
  metadata_filt_reorder<-filter_metadata(metadata,dat_filter_logscaled_rmnovar_t)
  
  fit_batch_regress=regress_batch(metadata_filt_reorder,dat_filter_logscaled_rmnovar_t,resids_file)
  zscore_batch_regress_all=get_z_scores(fit_batch_regress)
  myresids_batch_num=zscore_batch_regress_all[[1]];zscore_batch_regress=zscore_batch_regress_all[[2]]
  zscore_batch_regress_removedlowvar=remove_low_variance(metadata_filt_reorder,zscore_batch_regress,myresids_batch_num)
  print('alpha')
  reformat_and_write(zscore_batch_regress_removedlowvar,outfile)
  print('beta')
  # 
  # fit_sv_regress_all=regress_svs(metadata_filt_reorder,dat_filter_logscaled_rmnovar_t)
  # actual_svs=fit_sv_regress_all[[2]];fit_sv_regress=fit_sv_regress_all[[1]]
  # zscore_sv_regress_all=get_z_scores(fit_sv_regress)
  # myresids_sv_num=zscore_sv_regress_all[[1]];zscore_sv_regress=zscore_sv_regress_all[[2]]

  # fit_pc_regress_all=regress_pcs(metadata_filt_reorder,dat_filter_logscaled_rmnovar_t)
  # actual_pcs=fit_pc_regress_all[[2]];fit_pc_regress=fit_pc_regress_all[[1]]
  # zscore_pc_regress_all=get_z_scores(fit_pc_regress)
  # myresids_pc_num=zscore_pc_regress_all[[1]];zscore_pc_regress=zscore_pc_regress_all[[2]]
  # zscore_batch_regress_removedlowvar=remove_low_variance(metadata_filt_reorder,zscore_pc_regress,myresids_pc_num)
}

testing_things<-function(zscore_sv_regress,zscore_batch_regress){
  
  compare_zscores<-function(sample_id,gene){
    print(paste0("for sample ",sample_id, "/", gene, ' ----  ',"regressing batch: ",round(zscore_batch_regress[sample_id,gene],2), 
                 " ; sv only: ",round(zscore_sv_regress[sample_id,gene],2),
                 "; pc only: ",round(zscore_pc_regress[sample_id,gene],2)))
  }
  #BLOOD
  # compare_zscores('SAMPLE ID','GENE ID')

  
  library(ggplot2)
  #actual_svs=actual_svs%>%mutate(sample=rownames(acutal_svs))
  actual_svs_melted=melt(actual_svs$sv)
  colnames(actual_svs_melted)<-c("samples","svs","value")
  ggplot(actual_svs_melted,aes(y=samples,x=svs,fill=value))+
    geom_tile()+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    #geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient2(low="#43b0f0",mid = "#b8bcbf", high = "#de5414")+
    ggtitle('visualization of sva$sv')
  
  full_svs_mds=as.data.frame(modsv)%>%mutate(sample=rownames(modsv))
  colnames(full_svs_mds)<-gsub("RUNRUN",'',colnames(full_svs_mds))
  full_svs_mds_corr=cor(full_svs_mds[,2:16],full_svs_mds[,17:(ncol(full_svs_mds)-1)])
  full_svs_mds_corr_melted=melt(full_svs_mds_corr)
  ggplot(full_svs_mds_corr_melted,aes(x=Var2,y=Var1,fill=value))+
    geom_tile()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient2(low="#43b0f0",mid = "#b8bcbf", high = "#de5414",midpoint = 0)+
    ggtitle('Correlating SVs to Technical Md')
  full_svs_mds_corr_samples=cor(full_svs_mds[,18:ncol(full_svs_mds)-1],model.matrix(~sample,(full_svs_mds)))
  full_svs_mds_corr_samples=full_svs_mds_corr_samples%>%as.data.frame()%>%select_if(~any(abs(.) >= .2))%>%mutate(Var1=rownames(full_svs_mds_corr_samples))
  full_svs_mds_corr_samples_melted=melt(full_svs_mds_corr_samples) %>%mutate(variable=gsub('sample','',variable))%>%mutate(SV=gsub('SV','',Var1))
  ggplot(full_svs_mds_corr_samples_melted,aes(x=as.numeric(SV),y=variable,fill=value))+
    geom_tile()+
    xlab('SV')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient2(low="#43b0f0",mid = "#b8bcbf", high = "#de5414",midpoint = 0)+
    ggtitle('Correlating SVs to samples (FILTERED: abs(corr) > .2)')
  
  
  
  full_pcs_mds=as.data.frame(cbind(modsv[,2:16],actual_pcs$x[,1:30]))%>%mutate(sample=rownames(modsv))
  colnames(full_pcs_mds)<-gsub("RUNRUN",'',colnames(full_pcs_mds))
  full_pcs_mds_corr=cor(full_pcs_mds[,1:15],full_pcs_mds[,16:(ncol(full_pcs_mds)-1)])
  full_pcs_mds_corr_melted=melt(full_pcs_mds_corr)
  ggplot(full_pcs_mds_corr_melted,aes(x=Var2,y=Var1,fill=value))+
    geom_tile()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient2(low="#43b0f0",mid = "#b8bcbf", high = "#de5414",midpoint = 0)+
    ggtitle('Correlating pcs to Technical Md')
  full_pcs_mds_corr_samples=cor(full_pcs_mds[,17:ncol(full_pcs_mds)-1],model.matrix(~sample,(full_pcs_mds)))%>%as.data.frame()
  full_pcs_mds_corr_samples_filt=full_pcs_mds_corr_samples%>%select_if(~any((.)>0.5|.<(-.5),na.rm = T))%>%mutate(samples=rownames(full_pcs_mds_corr_samples))
  full_pcs_mds_corr_samples_melted=melt(full_pcs_mds_corr_samples_filt) %>%mutate(variable=gsub('sample','',variable))%>%mutate(PC=gsub('PC','',samples))
  ggplot(full_pcs_mds_corr_samples_melted,aes(x=as.numeric(PC),y=variable,fill=value))+
    geom_tile()+
    xlab('PC')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,  hjust=1))+
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient2(low="#43b0f0",mid = "#b8bcbf", high = "#de5414",midpoint = 0)+
    ggtitle('Correlating pcs to samples (FILTERED: abs(corr) > .2)')
  
  pc24=actual_pcs$rotation[,'PC24'][order(abs(actual_pcs$rotation[,'PC24']),decreasing = T)] 
  names(pc24)=gsub("\\..*",'',names(pc24))
  which(names(pc24)=='ENSG00000151327')
  pc12=actual_pcs$rotation[,'PC12'][order(abs(actual_pcs$rotation[,'PC12']),decreasing = T)] 
  names(pc12)=gsub("\\..*",'',names(pc12))
  which(names(pc12)=='ENSG00000140263')
  
  pc8=actual_pcs$rotation[,'PC8'][order(abs(actual_pcs$rotation[,'PC8']),decreasing = T)] 
  names(pc8)=gsub("\\..*",'',names(pc8))
  which(names(pc8)=='ENSG00000140263')
  pc10=actual_pcs$rotation[,'PC10'][order(abs(actual_pcs$rotation[,'PC10']),decreasing = T)] 
  names(pc10)=gsub("\\..*",'',names(pc10))
  which(names(pc10)=='ENSG00000140263')
  PC8plusPC10=abs(actual_pcs$rotation[,'PC8'])+abs(actual_pcs$rotation[,'PC10'])%>%as.data.frame()
  PC8plusPC10$gene=rownames(PC8plusPC10)
  colnames(PC8plusPC10)<-c('rotation','gene')
  PC8plusPC10=PC8plusPC10%>%mutate(rotation_rank=rank(-rotation))%>%mutate(gene=gsub("\\..*","",gene))
  PC8plusPC10%>%dplyr::filter(gene=='ENSG00000140263')
  # ggplot(,aes(=SV1,y=SV2,color=RIN))+geom_point()+scale_color_gradient2(low="#F7CA18",mid = "orange", high = "red",midpoint=7)+theme_bw()+
  #   geom_text(aes(label = round(value, 2)))
    

}

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
metadata_file<-as.character(opt$metadata)
tissue_type<-(as.character(opt$tissue))
outfile<-as.character(opt$outfile)
build<-as.character(opt$build)
outfile_transformed<-as.character(opt$outfile_transformed)
counts_file<-as.character(opt$counts)
num_svs<-as.character(opt$num_svs)
resids_file<-as.character(opt$resids_file)

main(metadata_file,tissue_type,counts_file,outfile,outfile_transformed,build,num_svs,resids_file)









