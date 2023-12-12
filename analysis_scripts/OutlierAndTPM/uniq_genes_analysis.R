library(data.table)
library(UpSetR)
library(ComplexHeatmap)
library(tidyverse)

###load needed data
md=fread("/Volumes/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv")
sample_to_status=md$affected_status
names(sample_to_status)=md$sample_id
sample_to_tiss=md$source_of_RNA
names(sample_to_tiss)=md$sample_id

uniq_genes_dir="/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/UniqGenesInfo"
hg38_uniq_genes=fread(paste0(uniq_genes_dir,"/s9_hg38_uniq_genes.csv"))
chm13_uniq_genes=fread(paste0(uniq_genes_dir,"/s11_extra_chm13_genes.csv"))

filter_tpms = function(df, people_filt = .2, tpm_filt = 0.15,min_sample_size=30) {
  genes = unique(df$gene)
  ind_filt_both = round(people_filt*length(unique(df$sample))) #this is 20% of people
  exp_df=(df)%>%group_by(gene)%>%mutate(n = sum(tpm >=tpm_filt))%>%dplyr::filter(n>ind_filt_both)
  return(exp_df%>%dplyr::select(-n))
}
###exp
read_in_file_hg38<-function(filename,uniq_genes){
  my_df=fread(filename)
  colnames(my_df)=c("sample","gene","transcript","tpm")
  my_df$gene<-gsub("\\..*","",my_df$gene)
  filter_df=filter_tpms(my_df)
  merged_df<-merge(filter_df,uniq_genes,by.x="gene",by.y="Gene ID")
  return(merged_df)
}
get_hg38_unique<-function(){
  
  hg38_fibroblast=read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
  hg38_blood=read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
  hg38_iPS=read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/iPS.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
  hg38_iPSC_NPC=read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/iPSC_NPC.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
  hg38_muscle=read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Muscle.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
  hg38_PBMC=read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/PBMC.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
  
  hg38_exp_allgenes<-rbind(
    cbind(hg38_fibroblast,tissue="fibroblast"),
    cbind(hg38_blood,tissue="blood"),
    cbind(hg38_iPS,tissue="iPS"),
    cbind(hg38_iPSC_NPC,tissue="iPSC_NPC"),
    cbind(hg38_muscle,tissue="muscle"),
    cbind(hg38_PBMC,tissue="PBMC"))
  
  
  
  # hg38_exp_genes_uniq=unique(rbind(unique(hg38_fibroblast[,c(1,5:9)]),
  #                                  unique(hg38_blood[,c(1,5:9)]),
  #                                  unique(hg38_iPS[,c(1,5:9)]),
  #                                  unique(hg38_iPSC_NPC[,c(1,5:9)]),
  #                                  unique(hg38_muscle[,c(1,5:9)]),
  #                                  unique(hg38_PBMC[,c(1,5:9)])))
  # hg38_exp_genes_uniq_num=rbind(fibroblast=nrow(unique(hg38_fibroblast[,c(1,5:9)])),
  #                               blood= nrow(unique(hg38_blood[,c(1,5:9)])),
  #                               iPS=nrow(unique(hg38_iPS[,c(1,5:9)])),
  #                               iPSC_NPC=nrow(unique(hg38_iPSC_NPC[,c(1,5:9)])),
  #                               muscle=nrow(unique(hg38_muscle[,c(1,5:9)])),
  #                               PBMC=nrow(unique(hg38_PBMC[,c(1,5:9)])))
  
  return(hg38_exp_allgenes)
  
}
hg38_exp_allgenes<-read.csv("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/hg38unique.genes.results")
ggplot(hg38_exp_genes_uniq,aes(x=`Gene biotype`))+geom_bar()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Genes unique to hg38")
ggplot(hg38_exp_genes_uniq,aes(x=`Unmapped Reason`))+geom_bar()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Genes unique to hg38")



length(unique(unlist(upset_list_hg38)))
upset_list_hg38<-(list(fibroblast=unique(hg38_fibroblast$gene),
     blood=unique(hg38_blood$gene),
     iPS=unique(hg38_iPS$gene),
     iPSC_NPC=unique(hg38_iPSC_NPC$gene),
     muscle=unique(hg38_muscle$gene),
     PBMC=unique(hg38_PBMC$gene)))
upset_matrix_hg38<-make_comb_mat(upset_list_hg38)
UpSet(upset_matrix_hg38)

read_in_file_chm13<-function(filename,uniq_genes){
  my_df=fread(filename)
  colnames(my_df)=c("sample","gene","transcript","tpm")
  my_df$gene<-gsub("\\..*","",my_df$gene)
  filter_df=filter_tpms(my_df)
  merged_df<-merge(filter_df,uniq_genes,by.x="gene",by.y="gene ID")
  return(merged_df)
  return(merged_df)
}

chm13_fibroblast=read_in_file_chm13("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13.dedupOptical_minMQ255_rsem.genes.results",chm13_uniq_genes)
chm13_blood=read_in_file_chm13("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13.dedupOptical_minMQ255_rsem.genes.results",chm13_uniq_genes)
chm13_iPS=read_in_file_chm13("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/iPS.chm13.dedupOptical_minMQ255_rsem.genes.results",chm13_uniq_genes)
chm13_iPSC_NPC=read_in_file_chm13("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/iPSC_NPC.chm13.dedupOptical_minMQ255_rsem.genes.results",chm13_uniq_genes)
chm13_muscle=read_in_file_chm13("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Muscle.chm13.dedupOptical_minMQ255_rsem.genes.results",chm13_uniq_genes)
chm13_PBMC=read_in_file_chm13("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/PBMC.chm13.dedupOptical_minMQ255_rsem.genes.results",chm13_uniq_genes)

chm13_exp_allgenes<-read.csv("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/chm13_unique_genes.genes.results")
chm13_exp_allgenes<-rbind(
  cbind(chm13_fibroblast,tissue="fibroblast"),
  cbind(chm13_blood,tissue="blood"),
  cbind(chm13_iPS,tissue="iPS"),
  cbind(chm13_iPSC_NPC,tissue="iPSC_NPC"),
  cbind(chm13_muscle,tissue="muscle"),
  cbind(chm13_PBMC,tissue="PBMC"))
chm13_exp_genes_uniq=chm13_exp_allgenes%>%group_by(gene)
chm13_exp_genes_uniq=unique(rbind(unique(chm13_fibroblast[,c(1,5:11)]),
                                 unique(chm13_blood[,c(1,5:11)]),
                                 unique(chm13_iPS[,c(1,5:11)]),
                                 unique(chm13_iPSC_NPC[,c(1,5:11)]),
                                 unique(chm13_muscle[,c(1,5:11)]),
                                 unique(chm13_PBMC[,c(1,5:11)])))
chm13_exp_genes_uniq_wtissue=unique(rbind(
                                 cbind(unique(chm13_fibroblast[,c(1,5:11)]),tissue="fibroblast"),
                                 cbind(unique(chm13_blood[,c(1,5:11)]),tissue="blood"),
                                 cbind(unique(chm13_iPS[,c(1,5:11)]),tissue="iPS"),
                                 cbind(unique(chm13_iPSC_NPC[,c(1,5:11)]),tissue="iPSC_NPC"),
                                 cbind(unique(chm13_muscle[,c(1,5:11)]),tissue="muscle"),
                                 cbind(unique(chm13_PBMC[,c(1,5:11)]),tissue="PBMC")))
chm13_exp_genes_uniq_num=rbind(fibroblast=nrow(unique(chm13_fibroblast[,c(1,5:11)])),
                          blood= nrow(unique(chm13_blood[,c(1,5:11)])),
                          iPS=nrow(unique(chm13_iPS[,c(1,5:11)])),
                          iPSC_NPC=nrow(unique(chm13_iPSC_NPC[,c(1,5:11)])),
                          muscle=nrow(unique(chm13_muscle[,c(1,5:11)])),
                          PBMC=nrow(unique(chm13_PBMC[,c(1,5:11)])))
chm13_exp_genes_uniq=chm13_exp_genes_uniq_wtissue%>%dplyr::select(-tissue)%>%unique()
nrow(chm13_exp_genes_uniq)
chm13_exp_genes_uniq%>%dplyr::filter(`Novel region`==1)%>%nrow()
chm13_exp_genes_uniq%>%dplyr::filter(`closest GENCODE ID`!="None")%>%nrow()
chm13_exp_genes_uniq%>%dplyr::filter(`Intersects medically relevant name list`==1) #%>%nrow()
chm13_exp_genes_uniq%>%dplyr::filter(`GRCh38 issue`==1)%>%nrow()
chm13_exp_genes_uniq%>%dplyr::filter(`biotype`=="protein_coding" |`biotype`=="lncRNA"  )%>%nrow()
ggplot(chm13_exp_genes_uniq_wtissue,aes(x=tissue,fill=biotype))+geom_bar()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Genes unique to chm13")

upset_list_chm13<-(list(fibroblast=unique(chm13_fibroblast$gene),
                       blood=unique(chm13_blood$gene),
                       iPS=unique(chm13_iPS$gene),
                       iPSC_NPC=unique(chm13_iPSC_NPC$gene),
                       muscle=unique(chm13_muscle$gene),
                       PBMC=unique(chm13_PBMC$gene)))
upset_list_chm13<-list(fibroblast=chm13_exp_allgenes%>%dplyr::filter(tissue=="fibroblast")%>%pull(gene)%>%unique(),
                       blood=chm13_exp_allgenes%>%dplyr::filter(tissue=="blood")%>%pull(gene)%>%unique(),
                       iPS=chm13_exp_allgenes%>%dplyr::filter(tissue=="iPS")%>%pull(gene)%>%unique(),
                       iPSC_NPC=chm13_exp_allgenes%>%dplyr::filter(tissue=="iPSC_NPC")%>%pull(gene)%>%unique(),
                       muscle=chm13_exp_allgenes%>%dplyr::filter(tissue=="muscle")%>%pull(gene)%>%unique(),
                       PBMC=chm13_exp_allgenes%>%dplyr::filter(tissue=="PBMC")%>%pull(gene)%>%unique())
upset_matrix_chm13<-make_comb_mat(upset_list_chm13)
UpSet(upset_matrix_chm13)


###exp violin plots
exp_allgenes=rbind(cbind(chm13_exp_allgenes[,c(1,2,4,13)],build="chm13"),
                   cbind(hg38_exp_allgenes[,c(1,2,4,10)],build="hg38"))
exp_allgenes$logTPM<-log10(exp_allgenes$tpm)
exp_allgenes=exp_allgenes%>%mutate(ifelse(logTPM=="-Inf",0,logTPM))
ggplot(exp_allgenes,aes(x=tissue,alpha=build,y=logTPM,fill=tissue))+
  geom_hline(yintercept=0.15,color="#bda039")+
  geom_violin(scale="width")+
  scale_fill_manual(values=c("#27A795","#D90025","#0470B5","#BD71DC","#AACA2F","#EF7C18"))+
  scale_alpha_manual(values=c(0.4,1))+
  ggtitle("Expression of build unique genes")+
  #ylim(c(0,5000))+
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  # scale_y_continuous(trans="log10")+
  theme_bw()

### reasons plot
upset_list_chm13<-list(fibroblast=chm13_exp_allgenes%>%dplyr::filter(tissue=="fibroblast")%>%pull(gene)%>%unique(),
                       blood=chm13_exp_allgenes%>%dplyr::filter(tissue=="blood")%>%pull(gene)%>%unique(),
                       iPS=chm13_exp_allgenes%>%dplyr::filter(tissue=="iPS")%>%pull(gene)%>%unique(),
                       iPSC_NPC=chm13_exp_allgenes%>%dplyr::filter(tissue=="iPSC_NPC")%>%pull(gene)%>%unique(),
                       muscle=chm13_exp_allgenes%>%dplyr::filter(tissue=="muscle")%>%pull(gene)%>%unique(),
                       PBMC=chm13_exp_allgenes%>%dplyr::filter(tissue=="PBMC")%>%pull(gene)%>%unique())

chm13_exp_genes_uniq<-chm13_exp_allgenes%>%dplyr::select(-tissue,-transcript,-tpm,-sample)%>%unique()
upset_reason_list_chm13<-(list(
  `non-syntenic novel region`=chm13_exp_genes_uniq%>%dplyr::filter(`Novel.region`==1)%>%pull(gene),
  `known GRCh38 issue`=chm13_exp_genes_uniq%>%dplyr::filter(`GRCh38.issue`==1)%>%pull(gene),
  `putative paralog`=chm13_exp_genes_uniq%>%dplyr::filter(`closest.GENCODE.ID`!='None')%>%pull(gene),
  `putative novel/unannotated gene`=chm13_exp_genes_uniq%>%dplyr::filter(biotype=='StringTie')%>%pull(gene),
  `no reason given`=chm13_exp_genes_uniq%>%dplyr::filter(biotype!='StringTie'&`closest GENCODE ID`=='None'&`GRCh38 issue`!=1&`Novel region`!=1)%>%pull(gene)
))
upset_reason_matrix_chm13<-make_comb_mat(upset_reason_list_chm13)
UpSet(upset_reason_matrix_chm13)

upset_reason_list_hg38<-(list(
  `chm13 duplication collapsed in hg38`=hg38_exp_genes_uniq%>%dplyr::filter(`Unmapped Reason`=="collapse")%>%pull(gene),
  `known GRCh38 issue`=chm13_exp_genes_uniq%>%dplyr::filter(`GRCh38 issue`==1)%>%pull(gene),
  `putative paralog`=chm13_exp_genes_uniq%>%dplyr::filter(`closest GENCODE ID`!='None')%>%pull(gene),
  `putative novel/unannotated gene`=chm13_exp_genes_uniq%>%dplyr::filter(biotype=='StringTie')%>%pull(gene),
  `no reason given`=chm13_exp_genes_uniq%>%dplyr::filter(biotype!='StringTie'&`closest GENCODE ID`=='None'&`GRCh38 issue`!=1&`Novel region`!=1)%>%pull(gene)
))
upset_reason_matrix_chm13<-make_comb_mat(upset_reason_list_chm13)
UpSet(upset_reason_matrix_chm13)



reason_hg38=hg38_exp_genes_uniq[,c("gene","Unmapped Reason")]%>%rename(reason=`Unmapped Reason`)%>%
  mutate(reason=ifelse(reason=="collapse","chm13 duplication collapsed in hg38", #low chm13 copy number
                       ifelse(reason=="no_alignment","no alignment to chm13",reason)
  ))
reason_chm13=chm13_exp_genes_uniq%>%
  mutate(neither_novel_nor_issue=`Novel region`+`GRCh38 issue`)%>%
  mutate(reason=ifelse(neither_novel_nor_issue==0,"neither_novel_nor_issue",
                       ifelse(`Novel region`==1, "novel_region","GRCh38_issue")))%>%
  dplyr::select(gene,reason)
all_reasons=rbind(cbind(reason_hg38,build="GRCh38"),
                  cbind(chm13_exp_genes_uniq[,c('gene','reason')],build="CHM13")) %>%
  mutate(build=factor(build,levels=c('GRCh38','CHM13')))
  # mutate(reason = factor(reason, 
  #                             levels=c("chm13 duplication collapsed in hg38","no alignment to chm13","hg38_error",
  #                               "novel_region","GRCh38_issue","neither_novel_nor_issue")))
ggplot(reason_hg38,aes(y=reason))+
  scale_fill_manual(values=c("#dea700","#01427a"))+
  # facet_wrap(~build,ncol = 1,scales = "free")+
  geom_bar(na.rm = T)+
  theme_bw()

###outliers
fibroblast_chm13eOutlier=fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%
  group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])
blood_chm13eOutlier=fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%
  group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])
fibroblast_hg38eOutlier=fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%
  group_by(sample_id)%>%dplyr::mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])
blood_hg38eOutlier=fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%
  group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])
hg38sOutlier=fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg38.sorted_dedupOptical_minMQ255.bed.gz") #%>%
colnames(hg38sOutlier)<-c("chr","start","end","cluster","sample_id","z",
                          "junc_name","score","strand","splice_site","acceptors_skipped","exons_skipped",
                          "donors_skipped","anchor","known_donor","known_junction","transcript","gene","x")
hg38sOutlier=hg38sOutlier%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(tissue=sample_to_tiss[sample_id])
blood_hg38sOutlier<-hg38sOutlier%>%dplyr::filter(tissue=="Blood")
fibroblast_hg38sOutlier<-hg38sOutlier%>%dplyr::filter(tissue=="Fibroblast")
chm13sOutlier=fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.chm13.sorted_dedupOptical_minMQ255.bed.gz") #%>%
colnames(chm13sOutlier)=c("chr","start","end","cluster","sample_id","z",
                  "transcript","gene","chm13_gene","collapsedon","x")
chm13sOutlier=chm13sOutlier%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(tissue=sample_to_tiss[sample_id])
blood_chm13sOutlier<-chm13sOutlier%>%dplyr::filter(tissue=="Blood")
fibroblast_chm13sOutlier<-chm13sOutlier%>%dplyr::filter(tissue=="Fibroblast")

##outlier plots
#####combine for rank plot
#get_outliers_in_atleast1<-function(fibroblast_chm13eOutlier,chm13_uniq_genes,by.y,uniq_build="chm13",tissue="fibroblast",exp_spl="exp"){
  outlier_genes=fibroblast_chm13eOutlier%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()
  fibroblast_eoutlier_chm13uniq<-cbind(merge(,chm13_uniq_genes,by.x="gene",by.y=by.y),uniq_build="chm13",tissue="fibroblast",exp_spl="exp")
  

fibroblast_eoutlier_chm13uniq<-cbind(merge(fibroblast_chm13eOutlier,chm13_uniq_genes,by.x="gene",by.y="gene ID"),uniq_build="chm13",tissue="fibroblast",exp_spl="exp")
blood_eoutlier_chm13uniq<-cbind(merge(blood_chm13eOutlier,chm13_uniq_genes,by.x="gene",by.y="gene ID"),uniq_build="chm13",tissue="blood",exp_spl="exp")
fibroblast_soutlier_chm13uniq<-cbind(merge(fibroblast_chm13sOutlier,chm13_uniq_genes,by.x="gene",by.y="gene ID"),uniq_build="chm13",tissue="fibroblast",exp_spl="spl")
blood_soutlier_chm13uniq<-cbind(merge(blood_chm13sOutlier,chm13_uniq_genes,by.x="gene",by.y="gene ID"),uniq_build="chm13",tissue="blood",exp_spl="spl")

fibroblast_eoutlier_hg38uniq<-cbind(merge(fibroblast_hg38eOutlier%>%mutate(gene=gsub("\\..*","",gene)),hg38_uniq_genes,by.x="gene",by.y="Gene ID"),uniq_build="hg38",tissue="fibroblast",exp_spl="exp")
blood_eoutlier_hg38uniq<-cbind(merge(blood_hg38eOutlier%>%mutate(gene=gsub("\\..*","",gene)),hg38_uniq_genes,by.x="gene",by.y="Gene ID"),uniq_build="hg38",tissue="blood",exp_spl="exp")
fibroblast_soutlier_hg38uniq<-cbind(merge(fibroblast_hg38sOutlier%>%mutate(gene=gsub("\\..*","",gene)),hg38_uniq_genes,by.x="gene",by.y="Gene ID"),uniq_build="hg38",tissue="fibroblast",exp_spl="spl")
blood_soutlier_hg38uniq<-cbind(merge(blood_hg38sOutlier%>%mutate(gene=gsub("\\..*","",gene)),hg38_uniq_genes,by.x="gene",by.y="Gene ID"),uniq_build="hg38",tissue="blood",exp_spl="spl")


uniq_eoutlier=rbind(blood_eoutlier_hg38uniq,fibroblast_eoutlier_hg38uniq)%>%dplyr::filter(abs(zscore)>3)

#rbind(blood_eoutlier_chm13uniq,fibroblast_eoutlier_chm13uniq)%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()%>%length()
#rbind(blood_eoutlier_hg38uniq,fibroblast_eoutlier_hg38uniq)%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()%>%length()
rbind(blood_soutlier_chm13uniq,fibroblast_soutlier_chm13uniq)[,-c(24)]%>%dplyr::filter(abs(z)>3)%>%pull(gene)%>%unique()%>%length()
rbind(blood_soutlier_hg38uniq,fibroblast_soutlier_hg38uniq)[,-c(24)]%>%dplyr::filter(abs(z)>3)%>%pull(gene)%>%unique()%>%length()

myranks=rbind(fibroblast_eoutlier_hg38uniq[,c(1,2,3,7,8,14,15,16)],
              blood_eoutlier_hg38uniq[,c(1,2,3,7,8,14,15,16)],
              fibroblast_soutlier_hg38uniq[,c(1,6,7,20,21,28,29,30)]%>%rename(zscore=z),
              blood_soutlier_hg38uniq[,c(1,6,7,20,21,28,29,30)]%>%rename(zscore=z),
              fibroblast_eoutlier_chm13uniq[,c(1,2,3,7,8,17,18,19)],
              blood_eoutlier_chm13uniq[,c(1,2,3,7,8,17,18,19)],
              fibroblast_soutlier_chm13uniq[,c(1,6,7,12,13,23,24,25)]%>%rename(zscore=z),
              blood_soutlier_chm13uniq[,c(1,6,7,12,13,23,24,25)]%>%rename(zscore=z))
ranks_test=myranks%>%dplyr::filter(rank_ind<20)
#ranks_test%>%dplyr::filter(status=="Case")%>%group_by(exp_spl,uniq_build,gene)%>%summarise(num_outliers=n())%>%group_by(exp_spl,uniq_build)%>%summarise(num_outliers=n())
ggplot(myranks,aes(x=rank_ind,fill=exp_spl))+
  facet_wrap(~uniq_build)+
  geom_bar() +
  scale_fill_manual(values=c("#88ace3","#4d7abf"))+
  xlab("Rank")+
  theme_classic()

num_outliers_per_sample=myranks%>%group_by(sample_id,uniq_build,tissue,exp_spl,status)%>%dplyr::filter(abs(zscore)>3)%>%summarise(num_outliers=n())
uniqoutliers_descriptive_plot=ggplot(num_outliers_per_sample%>%dplyr::filter(uniq_build!="chm13ensembl"),
                                 aes(x=exp_spl,y=num_outliers,fill=uniq_build))+
  geom_boxplot()+
  scale_alpha_manual(values=c(0.7,0.3))+
  scale_y_continuous(trans='log10')+
  #facet_wrap(~tissue,scale="free_y")+
  #scale_fill_manual(values=c("#dea700","#01427a"))+
  scale_fill_manual(values=c("#dea700","#01427a"))+
  ggtitle("Number of outliers per sample")+
  xlab("")+ylab("number of outliers")+
  theme_bw()
uniqoutliers_descriptive_plot

#####OR plot
get_OR<-function(fibroblast_hg38eOutlier,this_splexp,this_tissue){
  if(this_splexp=="expression"){
    mytable=fibroblast_hg38eOutlier%>%mutate(is_outlier=ifelse(abs(zscore)>3,"outlier","nonoutlier"))%>%
      ungroup()%>%dplyr::select(status,is_outlier)%>%table()
  }else{
    mytable=fibroblast_hg38eOutlier%>%mutate(is_outlier=ifelse(abs(z)>3,"outlier","nonoutlier"))%>%
      ungroup()%>%dplyr::select(status,is_outlier)%>%table()
  }
  
  fishers_mytable=fisher.test(mytable[,c("outlier","nonoutlier")])
  res=t(as.data.frame(c(pval=fishers_mytable$p.value,
                        odds_ratio=fishers_mytable$estimate,
                        lower_bound=fishers_mytable$conf.int[1],
                        upper_bound=fishers_mytable$conf.int[2])))
  rownames(res)=""
  return(cbind.data.frame(res,
                          tissue=this_tissue,
                          spl_or_exp=this_splexp))
}


##CHM13
get_outlier_info_chm13<-function(fibroblast_chm13eOutlier,chm13_uniq_genes,this_splexp,this_tissue){
  fibroblast_eoutlier_uniqgenes<-merge(fibroblast_chm13eOutlier,chm13_uniq_genes,by.x="gene",by.y="gene ID")
  outliers_here=fibroblast_eoutlier_uniqgenes%>%dplyr::filter(abs(zscore)>3)
  print(paste0("There are ",length(unique(outliers_here$gene))," unique outlier genes and ",
               length(outliers_here$gene), " outliers total"))
  outliers_here_info=unique(outliers_here[,c(1,9:ncol(outliers_here))])
  med_relevance= outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1)%>%nrow()
  novel_region= outliers_here_info%>%dplyr::filter(`Novel region`==1)%>%nrow()
  grch38_issue= outliers_here_info%>%dplyr::filter(`GRCh38 issue`==1)%>%nrow()
  novel_or_issue= outliers_here_info%>%dplyr::filter(`Either novel or issue`==1)%>%nrow()
  print(paste0("Number uniq outlier genes intersecting with..."))
  print(paste0("medically relevant name list: ", med_relevance))
  if(med_relevance>0){print(outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1))}
  print(paste0("novel or issue: ", novel_or_issue))
  print(paste0("novel region: ", novel_region))
  print(paste0("grch38 issue: ", grch38_issue))
  print(paste0("biotype: ")) 
  print(table(outliers_here_info$biotype))
  
  print(paste0("abs(zscore)>3 those in top 10 list=",
               outliers_here%>%dplyr::filter(rank_ind<11)%>%nrow()))
         print(outliers_here%>%dplyr::filter(rank_ind<11)%>%pull(status)%>%table)
  
  print(paste0(
               " and top 20=",
               outliers_here%>%dplyr::filter(rank_ind<21)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<21)%>%pull(status)%>%table)
  
  res=get_OR(fibroblast_eoutlier_uniqgenes,this_splexp,this_tissue)
  return(res)
}
get_outlier_info_splicing_chm13<-function(chm13sOutlier,chm13_uniq_genes,this_splexp,this_tissue){
  fibroblast_eoutlier_uniqgenes<-merge(chm13sOutlier,chm13_uniq_genes,by.x="chm13_gene",by.y="gene ID")
  outliers_here=fibroblast_eoutlier_uniqgenes%>%dplyr::filter(abs(z)>3)
  print(paste0("There are ",length(unique(outliers_here$gene))," unique outlier genes and ",
               length(outliers_here$gene), " outliers total"))
  outliers_here_info=unique(cbind(outliers_here[,1],outliers_here[,9:ncol(outliers_here)])) # sorry this is gross dunno freaks out otherwise
  med_relevance= outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1)%>%nrow()
  novel_region= outliers_here_info%>%dplyr::filter(`Novel region`==1)%>%nrow()
  grch38_issue= outliers_here_info%>%dplyr::filter(`GRCh38 issue`==1)%>%nrow()
  novel_or_issue= outliers_here_info%>%dplyr::filter(`Either novel or issue`==1)%>%nrow()
  print(paste0("Number uniq outlier genes intersecting with..."))
  print(paste0("medically relevant name list: ", med_relevance))
  if(med_relevance>0){print(outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1))}
  print(paste0("novel or issue: ", novel_or_issue))
  print(paste0("novel region: ", novel_region))
  print(paste0("grch38 issue: ", grch38_issue))
  print(paste0("biotype: ")) 
  print(table(outliers_here_info$biotype))
  print(paste0("abs(zscore)>3 those in top 10 list=",
               outliers_here%>%dplyr::filter(rank_ind<11)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<11)%>%pull(status)%>%table)
  print(paste0(
    " and top 20=",
    outliers_here%>%dplyr::filter(rank_ind<21)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<21)%>%pull(status)%>%table)
  
  res=get_OR(fibroblast_eoutlier_uniqgenes,this_splexp,this_tissue)
  return(res)
}
eoutlier_fibroblast_chm13uniq=get_outlier_info_chm13(fibroblast_chm13eOutlier,chm13_uniq_genes,"expression","fibroblast")
eoutlier_blood_chm13uniq=get_outlier_info_chm13(blood_chm13eOutlier,chm13_uniq_genes,"expression","blood")
soutlier_fibroblast_chm13uniq=get_outlier_info_splicing_chm13(fibroblast_chm13sOutlier,chm13_uniq_genes,"splicing","fibroblast")
soutlier_blood_chm13uniq=get_outlier_info_splicing_chm13(blood_chm13sOutlier,chm13_uniq_genes,"splicing","blood")

#case example known medical gene FCGR3B splicing
# blood_chm13eOutlier%>%dplyr::filter(gene=="LOFF_G0001310" & zscore>2)
# chm13sOutlier%>%dplyr::filter(chm13_gene=="LOFF_G0000174" & z>2)
# chm13sOutlier%>%dplyr::filter(chm13_gene=="LOFF_G0001426" & z>2)

###hg38 unique
get_outlier_info_hg38<-function(fibroblast_hg38eOutlier,hg38_uniq_genes,this_splexp,this_tissue){
  fibroblast_hg38eOutlier$gene<-sapply(strsplit(fibroblast_hg38eOutlier$gene,"\\."),"[[",1)
  fibroblast_eoutlier_uniqgenes<-merge(fibroblast_hg38eOutlier,hg38_uniq_genes,by.x="gene",by.y="Gene ID")
  outliers_here=fibroblast_eoutlier_uniqgenes%>%dplyr::filter(abs(zscore)>3)
  print(paste0("There are ",length(unique(outliers_here$gene))," unique outlier genes and ",
               length(outliers_here$gene), " outliers total"))
  outliers_here_info=unique(outliers_here[,c(1,9:ncol(outliers_here))])
  med_relevance= outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1)%>%nrow()
  unmapped_reason= outliers_here_info%>%dplyr::filter(`Unmapped Reason`==1)%>%nrow()
  print(paste0("Number uniq outlier genes intersecting with..."))
  print(paste0("medically relevant name list: ", med_relevance))
  if(med_relevance>0){print(outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1))}
  print(paste0("unmapped reason: ", unmapped_reason))
  print(paste0("biotype: ")) 
  print(table(outliers_here_info$`Gene biotype`))
  print(paste0("abs(zscore)>3 those in top 10 list=",
               outliers_here%>%dplyr::filter(rank_ind<11)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<11)%>%pull(status)%>%table)
  print(paste0(
    " and top 20=",
    outliers_here%>%dplyr::filter(rank_ind<21)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<21)%>%pull(status)%>%table)
  
  res=get_OR(fibroblast_eoutlier_uniqgenes,this_splexp,this_tissue)
  return(res)
}
get_outlier_info_splicing_hg38<-function(hg38sOutlier,hg38_uniq_genes,this_splexp,this_tissue){
  hg38sOutlier$gene<-sapply(strsplit(hg38sOutlier$gene,"\\."),"[[",1)
  sOutlier_uniqgenes<-merge(hg38sOutlier,hg38_uniq_genes,by.x="gene",by.y="Gene ID")
  outliers_here=sOutlier_uniqgenes%>%dplyr::filter(abs(z)>3)
  print(paste0("There are ",length(unique(outliers_here$gene))," unique outlier genes and ",
               length(outliers_here$gene), " outliers total"))
  outliers_here_info=unique(cbind(outliers_here[,1],outliers_here[,9:ncol(outliers_here)])) # sorry this is gross dunno freaks out otherwise
  med_relevance= outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1)%>%nrow()
  unmapped_reason= outliers_here_info%>%dplyr::filter(`Unmapped Reason`==1)%>%nrow()
  print(paste0("Number uniq outlier genes intersecting with..."))
  print(paste0("medically relevant name list: ", med_relevance))
  if(med_relevance>0){print(outliers_here_info%>%dplyr::filter(`Intersects medically relevant name list`==1))}
  print(paste0("unmapped reason: ", unmapped_reason))
  print(paste0("biotype: ")) 
  print(table(outliers_here_info$`Gene biotype`))
  print(paste0("abs(zscore)>3 those in top 10 list=",
               outliers_here%>%dplyr::filter(rank_ind<11)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<11)%>%pull(status)%>%table)
  print(paste0(
    " and top 20=",
    outliers_here%>%dplyr::filter(rank_ind<21)%>%nrow()))
  print(outliers_here%>%dplyr::filter(rank_ind<21)%>%pull(status)%>%table)
  
  res=get_OR(sOutlier_uniqgenes,this_splexp,this_tissue)
  return(res)
}

eoutlier_fibroblast_hg38uniq=get_outlier_info_hg38(fibroblast_hg38eOutlier,hg38_uniq_genes,"expression","fibroblast")
eoutlier_blood_hg38uniq=get_outlier_info_hg38(blood_hg38eOutlier,hg38_uniq_genes,"expression","blood")
soutlier_fibroblast_hg38uniq=get_outlier_info_splicing_hg38(fibroblast_hg38sOutlier,hg38_uniq_genes,"splicing","fibroblast")
soutlier_blood_hg38uniq=get_outlier_info_splicing_hg38(blood_hg38sOutlier,hg38_uniq_genes,"splicing","blood")

OR_fibroblast_exp_chm13=get_OR(fibroblast_chm13eOutlier,"expression","fibroblast")
OR_blood_exp_chm13=get_OR(blood_chm13eOutlier,"expression","blood")
OR_fibroblast_spl_chm13=get_OR(fibroblast_chm13sOutlier,"splicing","fibroblast")
OR_blood_spl_chm13=get_OR(blood_chm13sOutlier,"splicing","blood")
OR_fibroblast_exp_hg38=get_OR(fibroblast_hg38eOutlier,"expression","fibroblast")
OR_blood_exp_hg38=get_OR(blood_hg38eOutlier,"expression","blood")
OR_fibroblast_spl_hg38=get_OR(fibroblast_hg38sOutlier,"splicing","fibroblast")
OR_blood_spl_hg38=get_OR(blood_hg38sOutlier,"splicing","blood")


##PLOT ORS
all_OR=rbind(cbind(rbind(eoutlier_fibroblast_chm13uniq,eoutlier_blood_chm13uniq,soutlier_fibroblast_chm13uniq,soutlier_blood_chm13uniq),build="chm13",uniq="build unique"),
      cbind(rbind(eoutlier_fibroblast_hg38uniq,eoutlier_blood_hg38uniq,soutlier_fibroblast_hg38uniq,soutlier_blood_hg38uniq),build="hg38",uniq="build unique"),
      cbind(rbind(OR_fibroblast_exp_chm13,OR_blood_exp_chm13,OR_fibroblast_spl_chm13,OR_blood_spl_chm13),build="chm13",uniq="not unique"),
      cbind(rbind(OR_fibroblast_exp_hg38,OR_blood_exp_hg38,OR_fibroblast_spl_hg38,OR_blood_spl_hg38),build="hg38",uniq="not unique"))

all_OR=rbind(cbind(rbind(eoutlier_blood_chm13uniq,soutlier_blood_chm13uniq),build="chm13",uniq="build unique"),
             cbind(rbind(eoutlier_blood_hg38uniq,soutlier_blood_hg38uniq),build="hg38",uniq="build unique"),
             cbind(rbind(OR_blood_exp_chm13,OR_blood_spl_chm13),build="chm13",uniq="not unique"),
             cbind(rbind(OR_blood_exp_hg38,OR_blood_spl_hg38),build="hg38",uniq="not unique"))

OR_build_uniq=rbind(cbind(rbind(eoutlier_fibroblast_chm13uniq,eoutlier_blood_chm13uniq,soutlier_fibroblast_chm13uniq,soutlier_blood_chm13uniq),build="chm13",uniq="build unique"),
                    cbind(rbind(eoutlier_fibroblast_hg38uniq,eoutlier_blood_hg38uniq,soutlier_blood_hg38uniq),build="hg38",uniq="build unique"))

ggplot(OR_build_uniq%>%dplyr::filter(tissue=="blood"),
       aes(x=spl_or_exp,y=`odds_ratio.odds ratio`,color=build))+
  geom_point(aes(size = 2),position=position_dodge(width=.8)) +
  scale_color_manual(values=c("#dea700","#01427a"))+
  geom_hline(yintercept=1,color="#888888")+
  theme_bw()+
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),width=0.4, position = position_dodge(width=.8))+
  # labs(title = 'Odds an outlier is in a case as compared to a control', x = '', y = 'OR')
  ggtitle("Odds a build-unique gene outlier is in a case as compared to a control")+xlab("")+ylab("OR")



#FRG1
FRG1_chm13_genes=c("CHM13_G0043626","CHM13_G0034693","CHM13_G0055936","CHM13_G0034700","LOFF_G0001769","CHM13_G0036432","LOFF_G0001736","CHM13_G0015651",
                 "LOFF_G0001639","LOFF_G0002852","LOFF_G0000730","CHM13_G0056023","LOFF_G0000580","CHM13_G0034716","LOFF_G0000878","CHM13_G0034724",
                 "LOFF_G0000552","LOFF_G0000886","CHM13_G0055909","LOFF_G0000705","LOFF_G0000759","LOFF_G0001708","LOFF_G0001867","LOFF_G0002871")
FRG1_chm13_exp=chm13_exp_allgenes%>%dplyr::filter(gene %in% FRG1_chm13_genes)
outliers_FRG1BP4=blood_eoutlier_chm13uniq%>%dplyr::filter(gene=="LOFF_G0001639"& abs(zscore)>3)
blood_eoutlier_hg38<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))
blood_eoutlier_hg38$gene<-gsub("\\..*","",blood_eoutlier_hg38$gene)
blood_eoutlier_hg38%>%dplyr::filter(gene=="ENSG00000149531" & sample_id %in% outliers_FRG1BP4$sample_id )
FRG1BP4=chm13_exp_allgenes%>%dplyr::filter(gene=="LOFF_G0001639")

all_read_in_file_hg38<-function(filename,uniq_genes){
  my_df=fread(filename)
  colnames(my_df)=c("sample","gene","transcript","tpm")
  my_df$gene<-gsub("\\..*","",my_df$gene)
  filter_df=filter_tpms(my_df)
  #merged_df<-merge(filter_df,uniq_genes,by.x="gene",by.y="Gene ID")
  return(filter_df)
}
all_hg38_fibroblast=all_read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
all_hg38_blood=all_read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
all_hg38_iPS=all_read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/iPS.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
all_hg38_iPSC_NPC=all_read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/iPSC_NPC.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
all_hg38_muscle=all_read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Muscle.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)
all_hg38_PBMC=all_read_in_file_hg38("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/PBMC.hg38.dedupOptical_minMQ255_rsem.genes.results",hg38_uniq_genes)

all_hg38_exp_allgenes<-rbind(
  cbind(all_hg38_fibroblast,tissue="fibroblast"),
  cbind(all_hg38_blood,tissue="blood"),
  cbind(all_hg38_iPS,tissue="iPS"),
  cbind(all_hg38_iPSC_NPC,tissue="iPSC_NPC"),
  cbind(all_hg38_muscle,tissue="muscle"),
  cbind(all_hg38_PBMC,tissue="PBMC"))


FRG1BP<-all_hg38_exp_allgenes%>%dplyr::filter(gene=="ENSG00000149531")
FRG_genes=rbind(as.data.frame(FRG1BP[,c("gene","sample","tpm","tissue")]),
                as.data.frame(FRG1BP4[,c("gene","sample","tpm","tissue")]))
ggplot(FRG_genes,aes(x=tissue,y=tpm,fill=gene))+
  geom_violin()+
  theme_bw()+
  ggtitle("FRG1BP4 and FRG1BP")
#

