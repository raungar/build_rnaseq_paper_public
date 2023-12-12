library(data.table)
library(stringi) #chosen for supposed speediness
library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(STRINGdb)

get_sorted_residuals<-function(tpm_diff){
  tpm_diff_ranked<-tpm_diff%>% 
    dplyr::filter(!is.na(med_hg38) & !is.na(med_hg38)) %>%
    dplyr::filter(med_hg38!=0 & med_hg19 !=0)%>%
    mutate(rank_med_hg19 =rank(med_hg19,na.last="keep")) %>%
    mutate(rank_med_hg38 =rank(med_hg38,na.last="keep")) %>% as.data.frame()
  rownames(tpm_diff_ranked)<-as.character(c(1:nrow(tpm_diff_ranked)))
  
  lm_med<-(lm(rank_med_hg38~rank_med_hg19,data=tpm_diff_ranked))
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


tpm_diff_blood_gene<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.sorted_dedup_filtered.tpm_genesdiff.txt.gz")
tpm_diff_fibroblast_gene<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.sorted_dedup_filtered.tpm_genesdiff.txt.gz")

# tpm_diff_blood_iso<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparison/TPMDifferencesPrimary/Blood.sorted_dedup_filtered.tpm_isoformsdiff.txt.gz")
# tpm_diff_fibroblast_iso<-fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparison/TPMDifferencesPrimary/Fibroblast.sorted_dedup_filtered.tpm_isoformsdiff.txt.gz")


# library(biomaRt)
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# res <- getBM(attributes = c('ensembl_transcript_id_version', 
#                             'ensembl_gene_id', 
#                             'external_transcript_name',
#                             'external_gene_name'),
#              filters = 'ensembl_transcript_id_version', 
#              values = tpm_diff_fibroblast_iso$gene_id,
#              mart = mart)


fibroblast_tpm_diff_sortbyresid<-get_sorted_residuals(tpm_diff_fibroblast_gene)
blood_tpm_diff_sortbyresid<-get_sorted_residuals(tpm_diff_blood_gene)
muscle_tpm_diff_sortbyresid<-get_sorted_residuals(tpm_diff_muscle_gene)
fibroblast_tpm_diff_sortbyresid_iso<-get_sorted_residuals(tpm_diff_fibroblast_iso)
blood_tpm_diff_sortbyresid_iso<-get_sorted_residuals(tpm_diff_blood_iso)
#head(fibroblast_tpm_diff_sortbyresid_min0.5[,c("gene_id","med_hg19","med_hg38","rank_med_hg19","rank_med_hg38","predicted","resid","resid_adj")])
#head(blood_tpm_diff_sortbyresid[,c("gene_id","med_hg19","med_hg38","rank_med_hg19","rank_med_hg38","predicted","resid","resid_adj")],20)
head(blood_tpm_diff_sortbyresid_iso[,c("gene_id","med_hg19","med_hg38","rank_med_hg19","rank_med_hg38","resid_adj")],20)

fpr venn diagram
a=tpm_diff_blood_gene%>% dplyr::filter(is.na(med_hg38)) %>% pull(med_hg19) %>% table()
sum(a[names(a)>0.5])

##filt
tpm_diff_fibroblast_min0.5=tpm_diff_fibroblast%>% dplyr::filter(med_hg19>0.5 | med_hg38>0.5)
fibroblast_tpm_diff_sortbyresid_min0.5<-get_sorted_residuals(tpm_diff_fibroblast_min0.5)
tpm_diff_fibroblast_min0.5_iso=tpm_diff_fibroblast_iso%>% dplyr::filter(med_hg19>0.5 | med_hg38>0.5)
fibroblast_tpm_diff_sortbyresid_min0.5_iso<-get_sorted_residuals(tpm_diff_fibroblast_min0.5_iso)

tpm_diff_blood_min0.5=tpm_diff_blood%>% dplyr::filter(med_hg19>0.5 | med_hg38>0.5)
tpm_diff_blood_min0.5_iso=tpm_diff_blood_iso%>% dplyr::filter(med_hg19>0.5 | med_hg38>0.5)
blood_tpm_diff_sortbyresid_min0.5<-get_sorted_residuals(tpm_diff_blood_min0.5)
blood_tpm_diff_sortbyresid_min0.5_iso<-get_sorted_residuals(tpm_diff_blood_min0.5_iso)

tpm_diff_muscle_min0.5=tpm_diff_muscle%>% dplyr::filter(med_hg19>0.5 | med_hg38>0.5)
muscle_tpm_diff_sortbyresid_min0.5<-get_sorted_residuals(tpm_diff_muscle_min0.5)


###top percent
top_1percent_blood<- blood_tpm_diff_sortbyresid[1:quantile(1:nrow(blood_tpm_diff_sortbyresid),0.05),]
top_1percent_blood_min0.5<- blood_tpm_diff_sortbyresid_min0.5[1:quantile(1:nrow(blood_tpm_diff_sortbyresid),0.05),]
top_1percent_fibroblast<- fibroblast_tpm_diff_sortbyresid[1:quantile(1:nrow(fibroblast_tpm_diff_sortbyresid),0.01),]
top_1percent_fibroblast_min0.5<- fibroblast_tpm_diff_sortbyresid_min0.5[1:quantile(1:nrow(fibroblast_tpm_diff_sortbyresid),0.01),]
top_1percent_fibroblast_min0.5_iso<- fibroblast_tpm_diff_sortbyresid_min0.5_iso[1:quantile(1:nrow(fibroblast_tpm_diff_sortbyresid_iso),0.01),]
top_1percent_muscle<- muscle_tpm_diff_sortbyresid[1:quantile(1:nrow(muscle_tpm_diff_sortbyresid),0.01),]
top_1percent_muscle_min0.5<- muscle_tpm_diff_sortbyresid_min0.5[1:quantile(1:nrow(muscle_tpm_diff_sortbyresid),0.01),]



get_gostres<-function(top_1percent_blood,name_for_file){
  gostres <- gost(query =top_1percent_blood$gene_id, 
                  organism = "hsapiens", ordered_query = TRUE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL,  
                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  print(gostres$result)
  gostresforgem <- gost(query =top_1percent_blood$gene_id, evcodes = TRUE, multi_query = FALSE, 
                        sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"))
  
  gem <- gostresforgem$result[,c("term_id", "term_name", "p_value", "intersection")]
  colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
  gem$FDR <- gem$p.Val
  gem$Phenotype = "+1"
  gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
  write.table(gem, file = paste0("/Users/rachelungar/Documents/MontgomeryLab/RareDisease/Data/gProfiler_gem_",name_for_file,".txt"), sep = "\t", quote = F, row.names = F)
  return(gostres)
}
fibroblast_gostres<-get_gostres(top_1percent_fibroblast,"fibroblast")
fibroblast_gostres0.5<-get_gostres(top_1percent_fibroblast_min0.5,"fibroblast0.5")
fibroblast_gostres0.5_iso<-get_gostres(top_1percent_fibroblast_min0.5_iso,"fibroblast0.5_iso")
blood_gostres<-get_gostres(top_1percent_blood,"blood")
blood_gostres0.5<-get_gostres(top_1percent_blood_min0.5,"blood0.5")
muscle_gostres<-get_gostres(top_1percent_muscle,"muscle")
muscle_gostres0.5<-get_gostres(top_1percent_muscle_min0.5,"muscle0.5")

get_uniq_genes<-function(mygost,name_for_file){
  this_uniq_genes=c()
  for(ensg_line in (mygost$result$intersection)){
    these_genes=unique(unlist(strsplit(ensg_line,",")))
    this_uniq_genes=unique(c(this_uniq_genes,these_genes))
  }
  write.table(this_uniq_genes, file = paste0("/Users/rachelungar/Documents/MontgomeryLab/RareDisease/Data/gostplotgenes_",name_for_file,".txt"), sep = "\t", quote = F, row.names = F)
  return(this_uniq_genes)
}
blood_uniq_genes<-get_uniq_genes(blood_gostres,"blood")
blood_uniq_genes0.5<-get_uniq_genes(blood_gostres0.5,"blood0.5")
fibroblast_uniq_genes<-get_uniq_genes(fibroblast_gostres,"fibroblast")
fibroblast_uniq_genes0.5<-get_uniq_genes(fibroblast_gostres0.5,"fibroblast0.5")
fibroblast_uniq_genes0.5_iso<-get_uniq_genes(fibroblast_gostres0.5_iso,"fibroblast0.5_iso")
muscle_uniq_genes<-get_uniq_genes(muscle_gostres,"muscle")
muscle_uniq_genes0.5<-get_uniq_genes(muscle_gostres0.5,"muscle0.5")



# gostplot(gostres, capped = F, interactive = TRUE)
gostres<-fibroblast_gostres0.5_iso
publish_gosttable(gostres, highlight_terms = gostres$result,
                  use_colors = TRUE,
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)
gostplot(gostres)



### COMPARISON PLOT
plot_tpm_blood<-tpm_diff_blood
plot_tpm_blood[is.na(plot_tpm_blood)] <- as.numeric(-1)
plot_tpm_fibroblast<-tpm_diff_fibroblast
plot_tpm_fibroblast[is.na(plot_tpm_fibroblast)] <- as.numeric(-1)
plot_tpm_fibroblast_iso<-tpm_diff_fibroblast_iso
plot_tpm_fibroblast_iso[is.na(plot_tpm_fibroblast_iso)] <- as.numeric(-1)
plot_tpm_muscle<-tpm_diff_muscle
plot_tpm_muscle[is.na(plot_tpm_muscle)] <- as.numeric(-1)

ggplot(blood_tpm_diff_sortbyresid_iso,aes(x=rank_med_hg19,rank_med_hg38,color=abs(resid_adj)))+
  scale_color_gradient(low="blue",high="red")+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+ggtitle("blood hg19 vs. hg38 median TPM ranks (isoform)")+
  geom_point()# +xlim(0,100)+ylim(0,1500)

# top_5percent_resid_quantile=quantile(abs(lm_med_fibroblast$residuals),probs=c(0.95))
# top_1percent_resid_quantile=quantile(abs(lm_med_fibroblast$residuals),probs=c(0.99))
# top5percent_fibroblast<-tpm_diff_fibroblast_sorted[abs(lm_med_fibroblast$residuals)>top_5percent_resid_quantile,]
# top1percent_fibroblast<-tpm_diff_fibroblast_sorted[abs(lm_med_fibroblast$residuals)>top_1percent_resid_quantile,]

###SOLVED
solved=c("")

solved_df<-fread("/Volumes/groups/smontgom/shared/UDN/Data/solvedcases_udn_metricsreport_2022-02-10.csv")
solved_mim<-unique(unlist(str_split(solved_df$`Gene MIM number`,"((?![0-9]).)"))) #regex captures non-numerics
solved_mim_filt=solved_mim[which(str_length(solved_mim)>2)]
mim2gene<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/mim2gene_filt.txt")
colnames(mim2gene)<-c("mim","type","entrex","hgnc","ensg")
mim2ensg<-mim2gene$ensg
names(mim2ensg)<-mim2gene$mim
solved_ensg<-mim2ensg[solved_mim_filt] ##cant do 3 of them

solved_blood<-blood_tpm_diff_sortbyresid %>% dplyr::filter(gene_id %in% solved_ensg) ##33 in
percentile_blood=ecdf(abs(blood_tpm_diff_sortbyresid$resid_adj))
solved_fibroblast<-fibroblast_tpm_diff_sortbyresid %>% dplyr::filter(gene_id %in% solved_ensg) #32 in
percentile_fibroblast=ecdf(abs(fibroblast_tpm_diff_sortbyresid$resid_adj))
solved_muscle<-muscle_tpm_diff_sortbyresid %>% dplyr::filter(gene_id %in% solved_ensg) %>% arrange(desc(abs(resid_adj)))#32 in
percentile_muscle=ecdf(abs(muscle_tpm_diff_sortbyresid$resid_adj))

hist(percentile_fibroblast((abs(solved_fibroblast$resid_adj))),
     main = "residual quantile solved fibroblast",
     breaks = nrow(solved_fibroblast))
hist(percentile_blood((abs(solved_blood$resid_adj))),
     main = "residual quantile solved blood",
     breaks = nrow(solved_blood))

solved_all<-rbind(cbind(solved_blood,"tissue"="blood"),
                  cbind(solved_fibroblast,"tissue"="fibroblast"),
                  cbind(solved_muscle,"tissue"="muscle"))
ggplot(solved_all,aes(x=rank_med_hg38,y=rank_med_hg19,color=abs(resid)))+
  geom_abline(intercept=0,slope=1,color="red")+
  scale_color_gradient(low="#5de3ab",high="#280a47")+
  ggtitle("SOLVED ONLY ")+theme_bw()+
  geom_point(aes(alpha=0.6))+facet_grid(~tissue,scales="free_x")

top_notinhg19_butinhg38=tpm_diff_fibroblast_ranked%>%
  dplyr::filter(is.na(rank_med_hg19)) %>%
  dplyr::filter(med_hg38>15)
# top_notinhg38_butinhg19=tpm_diff_fibroblast_ranked%>% dplyr::filter(is.na(rank_med_hg38)) %>% dplyr::filter(med_hg19>15)

#color=abs(resid_adj),alpha=0.5)



#dplyr::filter(tpm_diff_fibroblast_sorted,med_hg38!=0 & med_hg19 !=0)


# tpm_diff_blood_sorted_filt<-tpm_diff_blood_sorted%>%dplyr::filter(abs(tpm_diff)>1)
manual_blood=tpm_diff_blood_sorted_filt %>% dplyr::filter(abs(tpm_diff)>1000)
manual_fibroblast=tpm_diff_fibroblast_sorted_filt %>% dplyr::filter(abs(tpm_diff)>1000)

library(effsize)
cohens_d_fibroblast<-(tpm_diff_fibroblast$mean_hg38-tpm_diff_fibroblast$mean_hg19)/
  sqrt(((tpm_diff_fibroblast$mean_hg38**2)+(tpm_diff_fibroblast$mean_hg19**2))/2)

found_genes<-c("")
#ggplot(tpm_diff_fibroblast,aes(x=gtex_or_blood,y=abs(MedZ)))+geom_violin()+ggtitle("Distribution of Median Z of each gene grouped by UDN/GTEx")+theme_bw()
ggplot(tpm_diff_fibroblast,aes(x=adjusted_diff))+geom_density() + theme_bw()+ggtitle("Fibroblast Median TPM Diff (HG38-HG19)") #+xlim(c(-1,1))
ggplot(tpm_diff_blood,aes(x=adjusted_diff))+geom_density() + theme_bw()+ggtitle("Blood Median TPM Diff (HG38-HG19)") #+xlim(c(-1,1))

lm_med_blood<-(lm(med_hg38~med_hg19,data=tpm_diff_blood))


