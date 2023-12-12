library(data.table)
library(stringi) #chosen for supposed speediness
library(ggplot2)
library(stringr)
library(tidyverse)
library(gprofiler2)

z_diff_blood_hg19<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz")
z_diff_fibroblast_hg19<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz")
z_diff_blood<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg38.chm13ensembl.txt.gz")
z_diff_fibroblast<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg38.chm13ensembl.txt.gz")
colnames(z_diff_blood)[3:4]<-colnames(z_diff_fibroblast)[3:4]<-c("z_hg38","z_chm13ensembl")
md=fread("/Volumes/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv")
sample_to_status=md$affected_status
names(sample_to_status)=md$sample_id

z_diff_blood$z_diff<-z_diff_blood$z_chm13-z_diff_blood$z_hg38
z_diff_blood_hg19$z_diff<-z_diff_blood_hg19$z_hg38-z_diff_blood_hg19$z_hg19
z_diff_fibroblast$z_diff<-z_diff_fibroblast$z_chm13-z_diff_fibroblast$z_hg38
z_diff_fibroblast_hg19$z_diff<-z_diff_fibroblast_hg19$z_hg38-z_diff_fibroblast_hg19$z_hg19
z_diff_blood_allbuilds<-merge(z_diff_blood%>%mutate(gene=str_replace(gene,"\\.[0-9]+","")),z_diff_blood_hg19,by=c("sample_id","gene"),all=T)%>%select(-z_hg38.y)
z_diff_fibroblast_allbuilds<-merge(z_diff_fibroblast%>%mutate(gene=str_replace(gene,"\\.[0-9]+","")),z_diff_fibroblast_hg19,by=c("sample_id","gene"))%>%select(-z_hg38.y)
z_diff_blood_allbuilds_melted<-melt(z_diff_blood_allbuilds)%>%dplyr::filter(!(variable%in%c("z_diff.x","z_diff.y")))
z_diff_fibroblast_allbuilds_melted<-melt(z_diff_fibroblast_allbuilds)%>%dplyr::filter(!(variable%in%c("z_diff.x","z_diff.y")))


this_gene="ENSG0000"
plot1=z_diff_blood%>%dplyr::filter(str_detect(gene,this_gene))%>%mutate(gene=str_replace(gene,"\\.[0-9]+(\\-[0-9])*",""))
to_plot<-merge(plot1,z_diff_blood_hg19%>%dplyr::filter(gene==this_gene),by=c("sample_id","gene"),all=T)%>%select(-z_hg38.y)
plotme=melt(to_plot)%>%dplyr::filter(!(variable%in%c("z_diff.x","z_diff.y"))) %>% group_by(sample_id,gene)%>%mutate(n=n())
plotme$variable=factor(plotme$variable,levels=c("z_hg19","z_hg38.x","z_chm13ensembl"))
ggplot(plotme%>%unique(),
       aes(x=variable,y=value,alpha=0.6))+geom_point()+theme_bw()+#geom_segment()+
  geom_line(aes(group=interaction(sample_id)))+ggtitle(this_gene)
  

blood_multt2t_genes=z_diff_blood%>%mutate(gene2=str_replace(gene,"-[0-9]",""))%>%
  group_by(sample_id,gene2)%>%mutate(n=n())%>%dplyr::filter(n>1)
fibroblast_multt2t_genes=z_diff_fibroblast%>%mutate(gene2=str_replace(gene,"-[0-9]",""))%>%
                          group_by(sample_id,gene2)%>%mutate(n=n())%>%dplyr::filter(n>1)#%>%pull(gene2)%>%unique()

z_diff_fibroblast_multt2t=z_diff_fibroblast%>%dplyr::filter(str_detect(gene,fibroblast_multt2t_genes))%>%mutate(ensembl_id=str_replace(gene,"-[0-9]",""))
z_diff_blood_multt2t=z_diff_blood%>%dplyr::filter(str_detect(gene,blood_multt2t_genes))%>%mutate(ensembl_id=str_replace(gene,"-[0-9]",""))

z_diff_fibroblast$z_diff<-z_diff_fibroblast$z_hg38-z_diff_fibroblast$z_hg19
z_diff_fibroblast<-z_diff_fibroblast%>%group_by(sample_id)%>%mutate(hg19_rank=rank(-abs(z_hg19)))%>%mutate(hg38_rank=rank(-abs(z_hg38)))%>%mutate(rankdiff=hg38_rank-hg19_rank)
fibroblast_multt2t_omim_genes=fibroblast_multt2t_genes%>%mutate(gene3=str_replace(gene2,"\\.[0-9]+",""))%>%dplyr::filter(gene3%in%mim2gene$ensg)
blood_multt2t_omim_genes=blood_multt2t_genes%>%mutate(gene3=str_replace(gene2,"\\.[0-9]+",""))%>%dplyr::filter(gene3%in%mim2gene$ensg)

ggplot(blood_multt2t_omim_genes%>%dplyr::filter(gene3=="ENSG00000"),
       aes(x=interaction(sample_id,gene2),y=z_chm13ensembl,alpha=0.2,color=z_diff))+geom_point()+ 
  scale_color_viridis_c(option = "plasma")+theme_bw()

blood_zdiff2=z_diff_blood%>%dplyr::filter((abs(z_chm13ensembl)>=3 & abs(z_hg38)<3) | (abs(z_hg38)>=3 & abs(z_chm13ensembl)<3))%>%dplyr::filter(abs(z_diff)>2) %>%group_by(gene)%>%mutate(n=n()) %>%arrange(desc(n))
#blood_zdiff2=z_diff_blood_hg19%>%dplyr::filter((abs(z_hg38)>=3 & abs(z_hg19)<3) | (abs(z_hg19)>=3 & abs(z_hg38)<3))%>%dplyr::filter(abs(z_diff)>2) %>%group_by(gene)%>%mutate(n=n()) %>%arrange(desc(n))
fibroblast_zdiff2=z_diff_fibroblast%>%dplyr::filter((abs(z_chm13ensembl)>=3 & abs(z_hg38)<3) | (abs(z_hg38)>=3 & abs(z_chm13ensembl)<3))%>%dplyr::filter(abs(z_diff)>2) %>%group_by(gene)%>%mutate(n=n()) %>%arrange(desc(n))
# tpm_diff_blood<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.sorted_dedup_filtered.tpm_genesdiff.txt.gz")
# tpm_diff_fibroblast<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.sorted_dedup_filtered.tpm_genesdiff.txt.gz")
# 
z_diff_blood%>%dplyr::filter((abs(z_hg19)>=3 & abs(z_hg38)<3) | (abs(z_hg38)>=3 & abs(z_hg19)<3))%>%dplyr::filter(abs(z_diff)>2)%>%pull(gene)%>%unique()%>%length()
z_diff_fibroblast%>%dplyr::filter((abs(z_hg19)>=3 & abs(z_hg38)<3) | (abs(z_hg38)>=3 & abs(z_hg19)<3))%>%dplyr::filter(abs(z_diff)>2)%>%pull(gene)%>%unique()%>%length()

blood_outliers_per_sample<-z_diff_blood%>%group_by(sample_id)%>%summarise(n_outliers_hg38=sum(abs(z_hg38)>3,na.rm = T),
                                                                       n_outliers_chm13=sum(abs(z_chm13ensembl)>3,na.rm=T))%>%mutate(n_outliers_diff=n_outliers_chm13-n_outliers_hg38)
fibroblast_outliers_per_sample<-z_diff_fibroblast%>%group_by(sample_id)%>%summarise(n_outliers_hg38=sum(abs(z_hg38)>3,na.rm = T),
                                                                          n_outliers_chm13=sum(abs(z_chm13ensembl)>3,na.rm=T))%>%mutate(n_outliers_diff=n_outliers_chm13-n_outliers_hg38)
ggplot(blood_outliers_per_sample,aes(x=n_outliers_diff))+geom_histogram()+theme_bw()+ggtitle("difference in number of outliers btwn chm13-hg38 (blood)")

get_sorted_residuals_z<-function(z_diff){
  z_diff_ranked<-z_diff%>% 
    dplyr::filter(!is.na(z_hg38) & !is.na(z_hg38)) %>%
    dplyr::filter(z_hg38!=0 & z_hg19 !=0)%>%
    mutate(rank_z_hg19 =rank(z_hg19,na.last="keep")) %>%
    mutate(rank_z_hg38 =rank(z_hg38,na.last="keep")) %>% as.data.frame()
  rownames(z_diff_ranked)<-as.character(c(1:nrow(z_diff_ranked)))
  
  lm_med<-(lm(rank_z_hg38~rank_z_hg19,data=z_diff_ranked))
  print(summary(lm_med))
  # studres<-studres(lm_med)
  # resid_sorted<-sort((abs(lm_med$residuals)),dec=T)
  # resid_sorted<-sort((abs(studres)),dec=T)
  predicted_vals<-summary(lm_med)$coefficients[1,1]+
    (summary(lm_med)$coefficients[2,1]*lm_med$fitted.values)
  z_diff_wresid<-cbind(z_diff_ranked,"predicted"=predicted_vals,resid=lm_med$residuals)%>%
    mutate(resid_adj=resid/log(predicted))
  #mutate(resid_adj=resid/sqrt(predicted))
  z_diff_sortbyresid<-z_diff_wresid %>% arrange(-abs(resid_adj))
  return(z_diff_sortbyresid)
}
get_sorted_residuals_tpm<-function(tpm_diff){
  tpm_diff_ranked<-tpm_diff%>% 
    dplyr::filter(!is.na(med_hg38) & !is.na(med_hg19)) %>%
    dplyr::filter(med_hg38!=0 & med_hg19 !=0)%>%
    mutate(rank_tpm_hg19 =rank(med_hg19,na.last="keep")) %>%
    mutate(rank_tpm_hg38 =rank(med_hg38,na.last="keep")) %>% as.data.frame()
  rownames(tpm_diff_ranked)<-as.character(c(1:nrow(tpm_diff_ranked)))
  
  lm_med<-(lm(rank_tpm_hg38~rank_tpm_hg19,data=tpm_diff_ranked))
  print(summary(lm_med))
  # studres<-studres(lm_med)
  # resid_sorted<-sort((abs(lm_med$residuals)),dec=T)
  # resid_sorted<-sort((abs(studres)),dec=T)
  predicted_vals<-summary(lm_med)$coefficients[1,1]+
    (summary(lm_med)$coefficients[2,1]*lm_med$fitted.values)
  tpm_diff_wresid<-cbind(tpm_diff_ranked,"predicted"=predicted_vals,resid=lm_med$residuals)%>%
    mutate(resid_adj=resid/log(predicted))
  #mutate(resid_adj=resid/sqrt(predicted))
  tpm_diff_sortbyresid<-tpm_diff_wresid %>% arrange(-abs(resid_adj))
  return(tpm_diff_sortbyresid)
}

fibroblast_tpm_diff_sortbyresid<-get_sorted_residuals_tpm(tpm_diff_fibroblast)
blood_tpm_diff_sortbyresid<-get_sorted_residuals_tpm(tpm_diff_blood)
fibroblast_z_diff_sortbyresid<-get_sorted_residuals_z(z_diff_fibroblast) #r2=0.7246
blood_z_diff_sortbyresid<-get_sorted_residuals_z(z_diff_blood) #r2=0.6849
z_diff_blood_sortbydiff<-z_diff_blood%>%arrange(desc(abs(z_diff)))
#head(fibroblast_z_diff_sortbyresid_min0.5[,c("gene_id","z_hg19","z_hg38","rank_z_hg19","rank_z_hg38","predicted","resid","resid_adj")])
#head(blood_z_diff_sortbyresid[,c("gene_id","z_hg19","z_hg38","rank_z_hg19","rank_z_hg38","predicted","resid","resid_adj")],20)
head(fibroblast_z_diff_sortbyresid[,c("gene","z_hg19","z_hg38","rank_z_hg19","rank_z_hg38","predicted","resid","resid_adj")],20)
fibroblast_z_diff_sortbyresid$z_diff<-fibroblast_z_diff_sortbyresid$z_hg38-fibroblast_z_diff_sortbyresid$z_hg19
blood_z_diff_sortbyresid$z_diff<-blood_z_diff_sortbyresid$z_hg38-blood_z_diff_sortbyresid$z_hg19

#get freq
z_diff_blood%>%dplyr::filter(abs(z_diff)>2) %>%group_by(gene)%>%mutate(n=n()) %>%arrange(desc(n))


#outliers in both/outliers in hg19
get_outlier_prop<-function(z, build_for_denominator,build_to_compare){
  print(paste0("What proportion of ",build_to_compare, " is similar in ",build_for_denominator,"?"))
  num_outliers_consistent=round(z%>%dplyr::filter(get(build_for_denominator)>=3 & get(build_to_compare)>=3)%>%nrow()/z%>%dplyr::filter(get(build_for_denominator)>=3)%>%nrow()*100,3)
  print(paste0("percent of outliers that are consistent: ",num_outliers_consistent,"%"))
  num_outliers_notseenother=round(z%>%dplyr::filter(get(build_for_denominator)>=3 & is.na(get(build_to_compare)))%>%nrow()/z%>%dplyr::filter(get(build_for_denominator)>=3)%>%nrow()*100,3)
  print(paste0("percent not seen in the other build: ",num_outliers_notseenother ,"%"))
  print(z%>%dplyr::filter(get(build_for_denominator)>=3 & is.na(get(build_to_compare)))%>%nrow()/z%>%dplyr::filter(get(build_for_denominator)>=3)%>%nrow())
  for(i in c(1,2,3)){
    denom=z%>%dplyr::filter(get(build_for_denominator)>=3)%>%nrow()
    numerator=(z%>%dplyr::filter(get(build_for_denominator)>=3 & get(build_to_compare)<i)%>%nrow())
    frac=numerator/denom*100
    print(paste0("z<",i, " is ", round(frac,3), "%"))
    
  }
}
get_outlier_prop(z_diff_fibroblast,"z_hg19","z_hg38")
get_outlier_prop(z_diff_fibroblast,"z_chm13ensembl","z_hg38")
get_outlier_prop(z_diff_blood,"z_hg38","z_chm13ensembl")
get_outlier_prop(z_diff_blood,"z_chm13ensembl","z_hg38")
##hg19
get_outlier_prop(z_diff_fibroblast_hg19,"z_hg19","z_hg38")
get_outlier_prop(z_diff_fibroblast_hg19,"z_hg38","z_hg19")
get_outlier_prop(z_diff_blood_hg19,"z_hg19","z_hg38")
get_outlier_prop(z_diff_blood_hg19,"z_hg38","z_hg19")
ggplot(blood_z_diff_sortbyresid,aes(x=abs(z_diff)))+geom_histogram()+theme_bw()+ggtitle("abs(zdiff) distr blood")+  
  scale_y_continuous(trans="log10")
ggplot(blood_z_diff_sortbyresid,aes(x=abs(z_diff)))+geom_density()+theme_bw()+ggtitle("abs(zdiff) distr blood")
ggplot(fibroblast_z_diff_sortbyresid,aes(x=abs(z_diff)))+geom_histogram()+theme_bw()+ggtitle("abs(zdiff) distr fibroblast")+  
  scale_y_continuous(trans="log10")

ggplot(z_diff_fibroblast,aes(x=abs(z_diff)))+geom_density() + 
  theme_bw()+ggtitle("Fibroblast Z Diff (chm13-hg38)") 
ggplot(z_diff_blood%,aes(x=adjusted_diff))+geom_density() + 
  theme_bw()+ggtitle("Blood Median TPM Diff (HG38-HG19)") #+xlim(c(-1,1))

ggplot(z_diff_blood%>%dplyr::filter(gene=="ENSG0000"),aes(x=z_hg38))+geom_histogram() +
  theme_bw()+xlim(c(-4.1,4.1))+
  geom_vline(xintercept =z_diff_blood%>%dplyr::filter(gene=="ENSG0000"&sample_id=="XXX")%>%pull(z_hg38))+
  ggtitle(" HG38") #+xlim(c(-1,1))

outlier_hg38<-z_diff_fibroblast%>%dplyr::filter(abs(z_hg38)>=3)
outlier_hg19<-z_diff_fibroblast%>%dplyr::filter(abs(z_hg19)>=3)
ggplot(outlier_onlyone,aes(x=abs(z_hg38)))+geom_density() + theme_bw()+
  ggtitle("z_hg19 for hg38 outliers") +xlim(c(-3,3))

##plots
# z_diff_fibroblast_melt<-melt(fibroblast_z_diff_sortbyresid,c("sample_id","gene","z_diff","resid_adj"))
med_abs_residadj<-abs(median(blood_tpm_diff_sortbyresid$resid_adj))
ggplot(fibroblast_tpm_diff_sortbyresid,#%>%dplyr::filter(abs(z_diff)>0.3)
       aes(x=med_hg19,y=med_hg38,color=abs(resid_adj),alpha=0.7))+
  scale_color_gradient2(high="#0d0630",mid="#94bfa7",midpoint =med_abs_residadj)+
  geom_abline (yintercept = 0,color="red",alpha=0.8)+  
  geom_point()+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(trans="log10")+
  theme_bw()+ggtitle("fibroblast hg19 vs. hg38 med TPM")
##z
ggplot(z_diff_fibroblast%>%dplyr::filter(abs(z_diff)>0.1),
       aes(x=z_hg38,y=z_chm13ensembl,color=abs(z_diff),alpha=0.7))+
  scale_color_gradient2(high="#0d0630",mid="#94bfa7",midpoint =0)+
  geom_abline (color="red",alpha=0.8)+  
  geom_point()+
  theme_bw()+ggtitle("fibroblast zscore hg38 vs chm13")


blood_tpm_combined_hg19<-fread("/Volumes/groups/smontgom/shared/UDN/PreprocessingHG19/eOutliers/Blood_combined_hg19.genes.results")
colnames(blood_tpm_combined_hg19)<-c("sample_id","gene","enst","tpm")
blood_tpm_combined_hg19$ensg<-sapply(strsplit(blood_tpm_combined_hg19$gene,"\\."),"[[",1)
gene_tpm_hg19<-blood_tpm_combined_hg19%>%dplyr::filter(ensg=="ENSG0000")
blood_tpm_combined_hg38<-fread("/Volumes/groups/smontgom/shared/UDN/PreprocessingHG38/eOutliers/Blood_combined_hg38.genes.results")
colnames(blood_tpm_combined_hg38)<-c("sample_id","gene","enst","tpm")
blood_tpm_combined_hg38$ensg<-sapply(strsplit(blood_tpm_combined_hg38$gene,"\\."),"[[",1)
gehe_tpm_hg38<-blood_tpm_combined_hg38%>%dplyr::filter(ensg=="ENSG0000")
ggplot(gehe_tpm_hg38,aes(x=tpm))+geom_histogram()+theme_bw(base_size=20)+
  #scale_x_continuous(trans='log10')+
  ggtitle(" hg38")+
  geom_vline(xintercept=PARG_tpm_hg38%>%dplyr::filter(sample_id=="XXXX")%>%pull(tpm))

gene_z<-blood_z_diff_sortbyresid%>%dplyr::filter(gene=="ENSG0000")
ggplot(gene_tpm_hg19,aes(x=z_hg19,y=z_hg38))+theme_bw()+
  #scale_x_continuous(trans='log10')+
  ggtitle(" hg19")+
  geom_abline (yintercept = 0,color="gray",alpha=0.8)+xlim(c(-4,4))+ylim(c(-4,4))+
  geom_hline(yintercept = 0 ,color="gray",alpha=0.8)+
  geom_vline(xintercept = 0,color="gray",alpha=0.8)+
  geom_point()+geom_point(data=gene_tpm_hg19%>%dplyr::filter(sample_id=="xxxx"),size=5,color="#5de3ab")


### COMPARISON PLOT
plot_z_blood<-z_diff_blood
plot_z_blood[is.na(plot_z_blood)] <- as.numeric(-1)
plot_z_fibroblast<-z_diff_fibroblast
plot_z_fibroblast[is.na(plot_z_fibroblast)] <- as.numeric(-1)
plot_z_muscle<-z_diff_muscle
plot_z_muscle[is.na(plot_z_muscle)] <- as.numeric(-1)

top_1percent_fibroblast<- fibroblast_z_diff_sortbyresid[1:quantile(1:nrow(fibroblast_z_diff_sortbyresid),0.01),]

ggplot(top_1percent_fibroblast,aes(x=rank_z_hg19,rank_z_hg38,color=abs(resid_adj)))+
  scale_color_gradient(low="blue",high="red")+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+ggtitle("muscle hg19 vs. hg38 median TPM ranks")+
  geom_point()# +xlim(0,100)+ylim(0,1500)

# top_5percent_resid_quantile=quantile(abs(lm_med_fibroblast$residuals),probs=c(0.95))
# top_1percent_resid_quantile=quantile(abs(lm_med_fibroblast$residuals),probs=c(0.99))
# top5percent_fibroblast<-z_diff_fibroblast_sorted[abs(lm_med_fibroblast$residuals)>top_5percent_resid_quantile,]
# top1percent_fibroblast<-z_diff_fibroblast_sorted[abs(lm_med_fibroblast$residuals)>top_1percent_resid_quantile,]

###SOLVED
solved=c("")
solved<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases.txt")
solved_df<-fread("/Volumes/groups/smontgom/shared/UDN/Data/solvedcases_udn_metricsreport_2022-02-10.csv") #426 uniq
solved_mim<-unique(unlist(str_split(solved_df$`Gene MIM number`,"((?![0-9]).)"))) #regex captures non-numerics
solved_mim_filt=solved_mim[which(str_length(solved_mim)>2)]
mim2gene<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/mim2gene_filt.txt")
colnames(mim2gene)<-c("mim","type","entrex","hgnc","ensg")
mim2ensg<-mim2gene$ensg
names(mim2ensg)<-mim2gene$mim
hgnc2ensg<-mim2gene$ensg
names(hgnc2ensg)<-mim2gene$hgnc
solved_ensg<-unique(mim2ensg[solved_mim_filt]) 
##OMIM
omim_blood<-blood_zdiff2%>%mutate(ensg=str_replace(gene,"\\..*$","")) %>%dplyr::filter(ensg %in% mim2gene$ensg)%>%mutate(case_status=sample_to_status[sample_id])
omim_fibroblast<-fibroblast_zdiff2 %>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(ensg %in% mim2gene$ensg)%>%mutate(case_status=sample_to_status[sample_id])
omim_blood%>%dplyr::filter(abs(z_hg38)>3 )%>%pull(gene)%>%unique()%>%length()
omim_blood%>%dplyr::filter(abs(z_chm13ensembl)>3)%>%pull(gene)%>%unique()%>%length()
ggplot(omim_blood,aes(x=z_hg38,z_chm13ensembl,color=case_status,alpha=.5))+theme_bw()+geom_point()
##SOLVED

solved_blood<-blood_zdiff2 %>% dplyr::filter(gene %in% solved_ensg) ##308
# percentile_blood=ecdf(abs(blood_z_diff_sortbyresid$resid_adj))
solved_fibroblast<-fibroblast_zdiff2 %>% dplyr::filter(gene %in% solved_ensg) #32 in
# percentile_fibroblast=ecdf(abs(fibroblast_z_diff_sortbyresid$resid_adj))
solved_muscle<-muscle_z_diff_sortbyresid %>% dplyr::filter(gene %in% solved_ensg) %>% arrange(desc(abs(resid_adj)))#32 in
percentile_muscle=ecdf(abs(muscle_z_diff_sortbyresid$resid_adj))

diffexp_blood_all_hg38_chm13=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01C_LIMMA_DiffExp/Blood.dedupOptical_minMQ255.diff_expression_chm13ensembl_vs_hg38.limma_voom_dream.all.txt.gz")
diffexp_fibroblast_allhg38_chm13=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01C_LIMMA_DiffExp/Fibroblast.dedupOptical_minMQ255.diff_expression_chm13ensembl_vs_hg38.limma_voom_dream.all.txt.gz")
diffexp_solvedUDNwide_blood=diffexp_blood_all_hg38_chm13%>%dplyr::filter(abs(logFC)>1&adj.P.Val<0.01)%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(ensg%in%solved_ensg)
diffexp_solvedUDNwide_fibroblast=diffexp_fibroblast_allhg38_chm13%>%dplyr::filter(abs(logFC)>1&adj.P.Val<0.01)%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(ensg%in%solved_ensg)

hist(percentile_fibroblast((abs(solved_fibroblast$resid_adj))),
     main = "residual quantile solved fibroblast",
     breaks = nrow(solved_fibroblast))
hist(percentile_blood((abs(solved_blood$resid_adj))),
     main = "residual quantile solved blood",
     breaks = nrow(solved_blood))

solved_all<-rbind(cbind(solved_blood,"tissue"="blood"),
                  cbind(solved_fibroblast,"tissue"="fibroblast")) #,
                  #cbind(solved_muscle,"tissue"="muscle"))
ggplot(solved_all,aes(x=med_hg19,y=med_hg38,color=abs(resid),shape=tissue,size=3))+
  geom_abline(intercept=0,slope=1,color="red")+
  scale_color_gradient(low="#5de3ab",high="#280a47")+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  ggtitle("Diagnostic/Solved Gene Expression/Rank across all UDN")+theme_bw(base_size = 18)+#xlim(c(0,260))+ylim(c(0,260))+
  geom_point(aes(alpha=0.8)) #+facet_grid(~tissue,scales="free_x")

top_notinhg19_butinhg38=z_diff_fibroblast_ranked%>%
  dplyr::filter(is.na(rank_z_hg19)) %>%
  dplyr::filter(z_hg38>15)
# top_notinhg38_butinhg19=z_diff_fibroblast_ranked%>% dplyr::filter(is.na(rank_z_hg38)) %>% dplyr::filter(z_hg19>15)

#color=abs(resid_adj),alpha=0.5)


solved_stanford_transcriptomics<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_blood_fibroblast.txt") 
rdid_to_ensg<-solved_stanford_transcriptomics$ENSG
names(rdid_to_ensg)<-solved_stanford_transcriptomics$RDID
solved_stanford_transcriptomics_blood<-solved_stanford_transcriptomics%>%dplyr::filter(Tissue=="Blood")
solved_stanford_transcriptomics_fibroblast<-solved_stanford_transcriptomics%>%dplyr::filter(Tissue=="Fibroblast")
#dplyr::filter(z_diff_fibroblast_sorted,z_hg38!=0 & z_hg19 !=0)
zdiff_solved_stanford_blood<-z_diff_blood%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(sample_id %in% names(rdid_to_ensg) & rdid_to_ensg[sample_id]==ensg)
zdiff_solved_stanford_fibroblast<-z_diff_fibroblast%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(sample_id %in% names(rdid_to_ensg) & rdid_to_ensg[sample_id]==ensg)
zdiff_solved_stanford_blood_hg19<-z_diff_blood_hg19%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(sample_id %in% names(rdid_to_ensg) & rdid_to_ensg[sample_id]==ensg)
zdiff_solved_stanford_fibroblast_hg19<-z_diff_fibroblast_hg19%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(sample_id %in% names(rdid_to_ensg) & rdid_to_ensg[sample_id]==ensg)

zdiff_solved_stanford_blood_allbuilds<-merge(zdiff_solved_stanford_blood,zdiff_solved_stanford_blood_hg19,by=c("sample_id","ensg"))
zdiff_solved_stanford_fibroblast_allbuilds<-merge(zdiff_solved_stanford_fibroblast,zdiff_solved_stanford_fibroblast_hg19,by=c("sample_id","ensg"))

zdiff_solved_stanford_blood_tpm<-blood_tpm_diff_sortbyresid%>%dplyr::filter(gene_id%in%rdid_to_ensg)
zdiff_solved_stanford_fibroblast<-fibroblast_tpm_diff_sortbyresid%>%dplyr::filter(gene_id%in%rdid_to_ensg)
#z_diff_blood
zdiff_solved_all<-rbind(cbind(zdiff_solved_stanford_blood_allbuilds,"tissue"="blood"),
                  cbind(zdiff_solved_stanford_fibroblast_allbuilds,"tissue"="fibroblast")) #,
ggplot(zdiff_solved_stanford_blood_hg19,aes(x=z_hg19,y=z_hg38,color=abs(z_diff)))+
  geom_abline(intercept=0,slope=1,color="red")+
  scale_color_gradient(low="#5de3ab",high="#280a47")+
  ggtitle("Solved Stanford Z-Scores ")+theme_bw(base_size = 18)+xlim(c(-4,4))+ylim(c(-4,4))+
  geom_point(aes(alpha=0.8),size=3) #+facet_grid(~tissue,scales="free_x")

###RANK
ggplot(zdiff_solved_all,aes(x=hg19_rank,y=hg38_rank,color=abs(rankdiff)))+
  geom_abline(intercept=0,slope=1,color="red")+
  scale_color_gradient(low="#5de3ab",high="#280a47")+
  ggtitle("Solved Z-Scores ")+theme_bw(base_size = 18) + #xlim(c(-4,4))+ylim(c(-4,4))+
  geom_point(aes(alpha=0.8),size=3) #+facet_grid(~tissue,scales="free_x")
melted_zdiff_solved_all<-melt(zdiff_solved_all%>%dplyr::filter(z_hg38!=0 | z_hg19!=0)%>%select(-z_hg19,-z_hg38,-z_diff,-rankdiff))
colnames(melted_zdiff_solved_all)<-c("sample_id","gene","tissue","build","rank")
ggplot(melted_zdiff_solved_all,aes(x=build,y=rank))+geom_violin()+
  scale_color_gradient(low="#5de3ab",high="#280a47")+
  geom_line(aes(group=sample_id))+
  scale_y_continuous(trans='log10')+
  ggtitle("Solved Rank Expression")+theme_bw(base_size = 18) + #xlim(c(-4,4))+ylim(c(-4,4))+
  geom_point(aes(alpha=0.8),size=3) #+facet_grid(~tissue,scales="free_x")



####BOOTSTRAP FOR SIGNIFICANCE
library(boot)
set.seed(235729)
get_resid_adj <- function(data, i){
  return(data[i,"resid_adj"])
}
boostrap_discreps_blood <- boot(blood_tpm_diff_sortbyresid,get_resid_adj,R=10000)
boostrap_discreps_blood_quantiles=quantile(boostrap_discreps_blood$t,c(0.025,0.975))
blood_tpm_diff_sortbyresid$is_sig<-ifelse(blood_tpm_diff_sortbyresid$resid_adj>boostrap_discreps_blood_quantiles[2] | blood_tpm_diff_sortbyresid$resid_adj<boostrap_discreps_blood_quantiles[1],
                                                  "sig", "notsig")
# 2.5%     97.5% 
#   -6.502758  2.307687 
boostrap_discreps_fibr <- boot(fibroblast_tpm_diff_sortbyresid,get_resid_adj,R=10000)
boostrap_discreps_fibr_quantiles=quantile(boostrap_discreps_fibr$t,c(0.025,0.975))
fibroblast_tpm_diff_sortbyresid$is_sig<-ifelse(fibroblast_tpm_diff_sortbyresid$resid_adj>boostrap_discreps_blood_quantiles[2] | fibroblast_tpm_diff_sortbyresid$resid_adj<boostrap_discreps_blood_quantiles[1],
                                          "sig", "notsig")
#v     2.5%     97.5% 
# -4.722902  1.629175 
solved_blood$is_sig<-ifelse(solved_blood$resid_adj>boostrap_discreps_blood_quantiles[2] |
                                                       solved_blood$resid_adj<boostrap_discreps_blood_quantiles[1],
                                                  "sig", "notsig")
solved_blood%>%dplyr::filter(is_sig=="sig")%>%nrow()
blood_tpm_diff_sortbyresid_discrep$is_sig<-ifelse(blood_tpm_diff_sortbyresid_discrep$resid_adj>boostrap_discreps_blood_quantiles[2] | blood_tpm_diff_sortbyresid_discrep$resid_adj<boostrap_discreps_blood_quantiles[1],
                                                  "sig", "notsig")
fibroblast_tpm_diff_sortbyresid_discrep$is_sig<-ifelse(fibroblast_tpm_diff_sortbyresid_discrep$resid_adj>boostrap_discreps_fibr_quantiles[2] | fibroblast_tpm_diff_sortbyresid_discrep$resid_adj<boostrap_discreps_fibr_quantiles[1],
                                                  "sig","notsig")
blood_sig_genes=blood_tpm_diff_sortbyresid_discrep%>%dplyr::filter(is_sig=="sig") %>%pull(gene_id)
discreps_genes_blood<-discreps_genes%>%
  dplyr::filter(ensg %in%blood_sig_genes)
fibr_sig_genes=fibroblast_tpm_diff_sortbyresid_discrep%>%dplyr::filter(is_sig=="sig") %>%pull(gene_id)
discreps_genes_fibr<-discreps_genes%>%
  dplyr::filter(ensg %in%fibr_sig_genes)


library(effsize)
cohens_d_fibroblast<-(z_diff_fibroblast$mean_hg38-z_diff_fibroblast$mean_hg19)/
  sqrt(((z_diff_fibroblast$mean_hg38**2)+(z_diff_fibroblast$mean_hg19**2))/2)

found_genes<-c("")
#ggplot(z_diff_fibroblast,aes(x=gtex_or_blood,y=abs(MedZ)))+geom_violin()+ggtitle("Distribution of Median Z of each gene grouped by UDN/GTEx")+theme_bw()
med_blood<-(lm(z_hg38~z_hg19,data=z_diff_blood))




###compare to exp diff

diffexp_blood=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01C_LIMMA_DiffExp/Blood_diff_expression.limma_voom_dream.significant.txt")
diffexp_fibroblast=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01C_LIMMA_DiffExp/Fibroblast_diff_expression.limma_voom_dream.significant.txt")

diffexp_blood_all=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01C_LIMMA_DiffExp/Blood_diff_expression.limma_voom_dream.all.txt")
diffexp_fibroblast_all=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01C_LIMMA_DiffExp/Fibroblast_diff_expression.limma_voom_dream.all.txt")



table(unique(str_replace(blood_zdiff2$gene,"\\-[0-9]",""))%in%(diffexp_blood%>%dplyr::filter(abs(logFC)>1)%>%pull(gene)));
table(unique(str_replace(fibroblast_zdiff2$gene,"\\-[0-9]",""))%in%(diffexp_fibroblast%>%dplyr::filter(abs(logFC)>1)%>%pull(gene)))

blood_diffexp_outliersThatChanged=diffexp_blood_all%>%dplyr::filter(gene %in% unique(blood_zdiff2$gene))
colnames(blood_diffexp_outliersThatChanged)[6]<-"padj"
blood_diffexp_outliersThatChanged=blood_diffexp_outliersThatChanged%>%mutate("significance"=ifelse(padj<0.05,"significant","not significant"))
fibroblast_diffexp_outliersThatChanged=diffexp_fibroblast_all%>%dplyr::filter(gene %in% unique(fibroblast_zdiff2$gene))
colnames(fibroblast_diffexp_outliersThatChanged)[6]<-"padj"
fibroblast_diffexp_outliersThatChanged=fibroblast_diffexp_outliersThatChanged%>%mutate("significance"=ifelse(padj<0.05,"significant","not significant"))

fibroblast_diffexp_outliersThatChanged%>%pull(significance)%>%table()
fibroblast_diffexp_outliersThatChanged%>%pull(significance)%>%table()
ggplot(data=fibroblast_diffexp_outliersThatChanged,aes(x=abs(logFC),fill=significance))+
  geom_histogram(bins=nrow(fibroblast_diffexp_outliersThatChanged)*.35)+theme_classic()+ scale_fill_manual(values=c("#cf5d06","#d1af06"))+xlab("differential expression abs(logFC)")+geom_vline(xintercept = 1,color="#807b6f")



###discrepent genes
discreps_genes<-fread("/Users/rachelungar/Documents/MontgomeryLab/RareDisease/Data/discreps.s4.txt")
colnames(discreps_genes)<-c("gene","num_variants_concordant","num_variants_hg19only","num_variants_hg38only",
                            "p_hg19","p_hg38","known_omim_gene","sig_gwas","dgd","psuedogene","ukbb")
genemap<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/genemap2.nopound.txt")
colnames(genemap)<-c("chr","start","end","cyto","computed_cyto","mim","hgnc","name","approved_symbol",
                     "entrez","ensembl","comments","pheno","mouse_id")
# genemap_spl <- strsplit(genemap$hgnc, split = ", ")
# genemap_renamed<-data.frame(ensg = rep(genemap_renamed$ensg, sapply(genemap_spl, length)), hgnc = unlist(genemap_spl))
# x=data.frame(V1 = rep(genemap_renamed$ensg, sapply(genemap_spl, length)), V2 = unlist(genemap_spl))
library("biomaRt")
mart <- useEnsembl(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl")
mybm=getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
           mart = mart)
hgnc_to_ensg<-mybm$ensembl_gene_id
names(hgnc_to_ensg)<-mybm$hgnc_symbol
hgnc_to_ensg["CBSL"]<-hgnc_to_ensg["CBS"]
hgnc_to_ensg["SMIM11A"]<-hgnc_to_ensg["SMIM11"]
hgnc_to_ensg["SMIM11B"]<-hgnc_to_ensg["SMIM11"]
discreps_genes$ensg<-hgnc_to_ensg[discreps_genes$gene]
diffexp_blood_discrep<-diffexp_blood%>%dplyr::filter(gene%in%discreps_genes$ensg)
diffexp_fibroblast_discrep<-diffexp_fibroblast%>%dplyr::filter(gene%in%discreps_genes$ensg)
outlier_blood_discrep<-blood_zdiff2%>%dplyr::filter(gene%in%discreps_genes$ensg)
outlier_fibroblast_discrep<-fibroblast_zdiff2%>%dplyr::filter(gene%in%discreps_genes$ensg)


blood_genes_outliers=blood_zdiff2%>%pull(gene)%>%unique()
blood_genes_outliers_expoverlap=blood_diffexp_outliersThatChanged%>%dplyr::filter(significance=="significant")%>%pull(gene)
blood_genes_outliers_discrepoverlap=diffexp_blood_discrep%>%pull(gene)
blood_genes_outliers_notinexpoverlap=blood_genes_outliers[!(blood_genes_outliers%in%blood_genes_outliers_expoverlap)] #17/169 htat are not in the 149 overlap of diffexp
table(blood_genes_outliers_notinexpoverlap%in%blood_genes_outliers_discrepoverlap) # are they in the discrep
table(blood_genes_outliers_discrepoverlap%in%blood_genes_outliers_expoverlap)
table(blood_genes_outliers_expoverlap%in%blood_genes_outliers_discrepoverlap) #number of outliers that overlap with exp in 

fibroblast_genes_outliers=fibroblast_zdiff2%>%pull(gene)%>%unique()
fibroblast_genes_outliers_expoverlap=fibroblast_diffexp_outliersThatChanged%>%dplyr::filter(significance=="significant")%>%pull(gene)
fibroblast_genes_outliers_discrepoverlap=diffexp_fibroblast_discrep%>%pull(gene)
fibroblast_genes_outliers_notinexpoverlap=fibroblast_genes_outliers[!(fibroblast_genes_outliers%in%fibroblast_genes_outliers_expoverlap)] #17/169 htat are not in the 149 overlap of diffexp
table(fibroblast_genes_outliers_notinexpoverlap%in%fibroblast_genes_outliers_discrepoverlap) # are they in the discrep
table(fibroblast_genes_outliers_discrepoverlap%in%fibroblast_genes_outliers_expoverlap)
table(fibroblast_genes_outliers_expoverlap%in%fibroblast_genes_outliers_discrepoverlap)


###venns?
library(nVennR)
fibroblast_venn=plotVenn(list("outliers"=fibroblast_genes_outliers,"discrepant genes"=discreps_genes$ensg,"differentially expressed" =diffexp_fibroblast$gene),
                         nCycles = 2000, labelRegions=F, fontScale=2)
blood_venn=plotVenn(list("outliers"=blood_genes_outliers,"discrepant genes"=discreps_genes$ensg,"differentially expressed" =diffexp_fibroblast$gene),
                         nCycles = 2000, labelRegions=F, fontScale=2)
fibroblast_venn=plotVenn(list("outliers"=fibroblast_genes_outliers,"discrepant genes"=discreps_genes$ensg,
                              "differentially expressed" =diffexp_fibroblast_all%>%dplyr::filter(adj.P.Val<0.05&abs(logFC)>1)%>%pull(gene)%>%unique()),
                         nCycles = 2000, labelRegions=F, fontScale=2)
blood_venn=plotVenn(list("outliers"=blood_genes_outliers,"discrepant genes"=discreps_genes$ensg,
                              "differentially expressed" =diffexp_blood_all%>%dplyr::filter(adj.P.Val<0.05&abs(logFC)>1)%>%pull(gene)%>%unique()),
                    nCycles = 2000, labelRegions=F, fontScale=2)


###t2t uniq
s11genes=fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/CHM13v2.0/chm13.T2T.s11.extragenes.txt",header=F)%>%pull(V1)
no_ensemblgenes=fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/CHM13v2.0/chm13.T2T_unique_genes.txt")%>%pull(gene_id)
fibroblast_chm13eOutlier=fread("/Volumes/groups/smontgom/shared/UDN/PreprocessingT2T/eOutliers/Fibroblast.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))
blood_chm13eOutlier=fread("/Volumes/groups/smontgom/shared/UDN/PreprocessingT2T/eOutliers/Blood.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))

#num of these genes expressed
fibroblast_eoutlier_s11genes<-fibroblast_chm13eOutlier%>%dplyr::filter(gene%in%s11genes)
fibroblast_eoutlier_no_ensemblgenes<-fibroblast_chm13eOutlier%>%dplyr::filter(gene%in%no_ensemblgenes)
blood_eoutlier_s11genes<-blood_chm13eOutlier%>%dplyr::filter(gene%in%s11genes)
blood_eoutlier_no_ensemblgenes<-blood_chm13eOutlier%>%dplyr::filter(gene%in%no_ensemblgenes)
print(paste0(length(unique(fibroblast_eoutlier_s11genes$gene))," of ",length(unique(s11genes))," genes expressed in fibroblast"))
print(paste0(length(unique(blood_eoutlier_s11genes$gene))," of ",length(unique(s11genes))," genes expressed in blood"))


fibroblast_eoutlier_s11genes%>%dplyr::filter(abs(zscore)>3)%>%nrow(); fibroblast_eoutlier_s11genes%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()%>%length(); fibroblast_eoutlier_s11genes%>%dplyr::filter(abs(zscore)>3&rank_ind<11)%>%nrow()
fibroblast_eoutlier_no_ensemblgenes%>%dplyr::filter(abs(zscore)>3)%>%nrow();fibroblast_eoutlier_no_ensemblgenes%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()%>%length(); fibroblast_eoutlier_no_ensemblgenes%>%dplyr::filter(abs(zscore)>3&rank_ind<11)%>%nrow()
blood_eoutlier_s11genes%>%dplyr::filter(abs(zscore)>3)%>%nrow();blood_eoutlier_s11genes%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()%>%length(); blood_eoutlier_s11genes%>%dplyr::filter(abs(zscore)>3&rank_ind<11)%>%nrow()
blood_eoutlier_no_ensemblgenes%>%dplyr::filter(abs(zscore)>3)%>%nrow();blood_eoutlier_no_ensemblgenes%>%dplyr::filter(abs(zscore)>3)%>%pull(gene)%>%unique()%>%length(); blood_eoutlier_no_ensemblgenes%>%dplyr::filter(abs(zscore)>3&rank_ind<11)%>%nrow()

#top genes
blood_eoutlier_s11genes%>%dplyr::filter(rank_ind<=10)%>%mutate(status=sample_to_status[sample_id])%>%dplyr::filter(status=="Case")


##amount of genes being expressed
