library("biomaRt")
library(data.table)
library(tidyverse)
library(ggpubr)

md=fread("/Volumes/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv")
sample_to_status=md$affected_status
names(sample_to_status)=md$sample_id
sample_to_udnid=md$indv_id
names(sample_to_udnid)=md$sample_id
sample_to_tiss=md$source_of_RNA
names(sample_to_tiss)=md$sample_id
sample_to_status=md$affected_status
names(sample_to_status)=md$sample_id


chrY_chm13_genes=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/chm13_chrY.txt",header = F)
colnames(chrY_chm13_genes)<-c("hgnc","ensembl","chm13_id")
chrY_chm13_genes<-chrY_chm13_genes%>%mutate(gene=gsub("\\..*","",ensembl))



###Load everything in
eoutliers_blood_hg19<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Blood.hg19.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_fibroblast_hg19<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Fibroblast.hg19.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_blood_hg38<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_fibroblast_hg38<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_blood_chm13<-fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_fibroblast_chm13<-fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_blood_chm13ensembl<-fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13ensembl.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
eoutliers_fibroblast_chm13ensembl<-fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13ensembl.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"))%>%mutate(status=sample_to_status[sample_id])%>%mutate(ensg=str_replace(gene,"\\..*$",""))

soutliers_hg19<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg19/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg19.sorted_dedupOptical_minMQ255.bed.gz") #%>%mutate(status=sample_to_status["sample"])
soutliers_hg38<-fread("/Volumes/groups/smontgom/shared/UDN/Output/hg38/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg38.sorted_dedupOptical_minMQ255.bed.gz") #%>%mutate(status=sample_to_status["sample"])
soutliers_chm13<-fread("/Volumes/groups/smontgom/shared/UDN/Output/chm13/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.chm13.sorted_dedupOptical_minMQ255.bed.gz") #%>%mutate(status=sample_to_status["sample"])
colnames(soutliers_hg19)<-c("chr","pos1","pos2","clus","sample","z","enst","gene","x") #,"status")
colnames(soutliers_hg38)<-c("chr","pos1","pos2","clus","sample","z","junc","junc_support","strand","splicesite","a","b","c","noveltype","novel_donor","novel_acceptor","enst","gene","x") #,"status")
colnames(soutliers_chm13)<-c("chr","pos1","pos2","clus","sample","z","enst","gene","chm13_id","gene2","x") #,"status")
soutliers_hg19=soutliers_hg19%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
soutliers_hg38=soutliers_hg38%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
soutliers_chm13=soutliers_chm13%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])%>%mutate(ensg=str_replace(gene,"\\..*$",""))
soutliers_chm13ensembl=soutliers_chm13%>%dplyr::filter(str_detect(gene2,"ENSG"))
soutliers_fibroblast_hg19=soutliers_hg19%>%dplyr::filter(tissue=="Fibroblast")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_blood_hg19=soutliers_hg19%>%dplyr::filter(tissue=="Blood")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_fibroblast_hg38=soutliers_hg38%>%dplyr::filter(tissue=="Fibroblast")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_blood_hg38=soutliers_hg38%>%dplyr::filter(tissue=="Blood")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_fibroblast_chm13=soutliers_chm13%>%dplyr::filter(tissue=="Fibroblast")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_blood_chm13=soutliers_chm13%>%dplyr::filter(tissue=="Blood")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_fibroblast_chm13ensembl=soutliers_chm13ensembl%>%dplyr::filter(tissue=="Fibroblast")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))
soutliers_blood_chm13ensembl=soutliers_chm13ensembl%>%dplyr::filter(tissue=="Blood")%>%mutate(rank_ind=rank(desc(abs(z)),ties.method="first"))


##Figure A
get_outlier_summary_stat<-function(my_varname,sample_to_status){
  mydata=get(my_varname)
  tissue=sapply(strsplit(my_varname,"_"),"[[",2)
  build=sapply(strsplit(my_varname,"_"),"[[",3)
  if(grepl("^e",my_varname)){exp_spl="expression"}else{exp_spl="splicing"}
  if(exp_spl=="expression"){zname="zscore"; sample_name="sample_id"}else{zname="z";sample_name="sample"}
  if(exp_spl=="expression"){
    over=mydata%>%dplyr::filter(get(zname)>3)%>%group_by(sample_id)%>%summarise(num_outlier=n())%>%mutate(status=sample_to_status[sample_id])
    
      under=mydata%>%dplyr::filter(get(zname)<(-3))%>%group_by(sample_id)%>%summarise(num_outlier=n())%>%mutate(status=sample_to_status[sample_id])
      return(as.data.frame(rbind(cbind(exp_spl,tissue,build,over_under="over",over),
                                 cbind(exp_spl,tissue,build,over_under="under",under))))
  }else{
    over=mydata%>%dplyr::filter(get(zname)>3)%>%group_by(sample)%>%summarise(num_outlier=n())%>%mutate(status=sample_to_status[sample])%>%rename(sample_id=sample)
    
    return(as.data.frame(cbind(exp_spl,tissue,build,over_under="over",over)))
  }

}
# 
# myzscores=c("eoutliers_blood_hg19","eoutliers_blood_hg38","eoutliers_blood_chm13","eoutliers_blood_chm13ensembl",
#             "eoutliers_fibroblast_hg19","eoutliers_fibroblast_hg38","eoutliers_fibroblast_chm13","eoutliers_fibroblast_chm13ensembl",
            myzscores=c("soutliers_blood_hg19","soutliers_blood_hg38","soutliers_blood_chm13","soutliers_blood_chm13ensembl",
            "soutliers_fibroblast_hg19","soutliers_fibroblast_hg38","soutliers_fibroblast_chm13","soutliers_fibroblast_chm13ensembl")
summ_stats_outliers_df=data.frame()
for(this_zscore in myzscores){
  print(this_zscore)
  summ_stats_outliers_df=rbind(summ_stats_outliers_df,
                               get_outlier_summary_stat(this_zscore,sample_to_status))
}

outliers_descriptive_plot=ggplot(summ_stats_outliers_df%>%dplyr::filter(build!="chm13ensembl"),
                                 aes(x=interaction(exp_spl,over_under),y=num_outlier,fill=build))+
  geom_boxplot()+
  # scale_alpha_manual(values=c(0.7,0.3))+
  scale_y_continuous(trans='log10')+
  #facet_wrap(~tissue,scale="free_y")+
  scale_fill_manual(values=c("#ab1b1b","#dea700","#1064ad"))+
  # scale_color_manual(values=c("#8C1515","#dea700","#01427a"))+
  scale_color_manual(values=c("green"))+
  ggtitle("Number of outliers per sample")+
  xlab("")+ylab("number of outliers")+
  scale_x_discrete(labels=c("Under eOutlier", "Over eOutlier", "sOutlier"))+
  theme_bw()
outliers_descriptive_plot


###solved cases
get_omim<-function(min_confidence){
  mart <- useEnsembl(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl")
  mybm=getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             mart = mart,
             useCache = FALSE)
  hgnc_to_ensg<-mybm$ensembl_gene_id
  names(hgnc_to_ensg)<-mybm$hgnc_symbol
  ensg_to_hgnc<-mybm$hgnc_symbol
  names(ensg_to_hgnc)<-mybm$ensembl_gene_id
  ensg_to_hgnc=ensg_to_hgnc[(ensg_to_hgnc!="")]
  hgnc_more=as.data.frame(fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/hgnc_convert.txt"))
  ensg_to_hgnc_more=hgnc_more[,10]
  names(ensg_to_hgnc_more)=(hgnc_more[,2])
  combined_ensg_to_hgnc=unique(rbind(stack(ensg_to_hgnc_more),stack(hgnc_to_ensg)))
  colnames(combined_ensg_to_hgnc)<-c("ensg","hgnc")
  # write_tsv(combined_ensg_to_hgnc, "/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/hgnc_ensg_convert_complete.txt")
  ref_dir="/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/"
  pheno2gene="morbidmap.txt"
  pheno = fread(paste0(ref_dir, pheno2gene))
  colnames(pheno) = c("pheno","gene","mim","cyto")
  mim2gene=fread(paste0(ref_dir,"mim2gene.txt"))
  mim_to_ensg<-mim2gene$`Ensembl Gene ID (Ensembl)`
  names(mim_to_ensg)<-as.factor(mim2gene$`# MIM Number`)
  
  pheno_clean = pheno %>% 
    extract(pheno, into = c("pheno", "confidence"), "(.*) \\(([0-9]+)\\)$") %>%
    mutate(confidence = as.numeric(confidence)) %>% 
    extract(pheno, into = c("pheno", "pheno_mim"), "(.*), ([0-9]{6})$") %>% 
    mutate(status = ifelse(confidence == 1, "association", 
                           ifelse(confidence == 2, "linakge_mapping", 
                                  ifelse(confidence == 3, "known_molec_cause", 
                                         ifelse(confidence == 4, "del_dup_syndrome", NA))))) %>% 
    separate_rows(gene, convert = T, sep=", ") %>% 
    mutate(ensg = mim_to_ensg[mim]) %>% 
    mutate(ensg=ifelse(ensg=="" | is.na(ensg),hgnc_to_ensg[gene],ensg))%>%
    mutate(ensg=ifelse(ensg=="" | is.na(ensg),ensg_to_hgnc_more[gene],ensg))%>%
    
    distinct() #%>%  #
  # add something to fill down ensg within groups
  # group_by(pheno, pheno_mim, confidence, status, cyto) %>% 
  # mutate(approved_symbol = ifelse(!is.na(ensg), gene, NA)) %>% 
  # summarize(ensg = toString(na.omit(ensg)),
  #           approved_symbol = toString(na.omit(approved_symbol)),
  #           hgnc = toString(na.omit(gene))) %>% 
  # ungroup() %>% 
  # separate_rows(ensg, approved_symbol, convert=T, sep=", ") %>% 
  # group_by(ensg, hgnc, approved_symbol) %>% 
  # summarize(pheno = toString(unique(na.omit(pheno))),
  #           pheno_mim  = toString(unique(na.omit(pheno_mim ))),
  #           confidence = toString(unique(na.omit(confidence))),
  #           status = toString(unique(na.omit(status)))) %>% 
  # ungroup()
  pheno_filtered<-pheno_clean%>%dplyr::filter(confidence>=min_confidence)
  return(pheno_filtered)
}
omim_df=get_omim(min_confidence=3)
omim_genes=unique(omim_df$ensg)
omim_dic=as.factor(omim_df$ensg)
names(omim_dic)=omim_df$gene
solved_and_candidates<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_AND_candidate_dec2022.txt")
get_solved_samples<-function(myvar,solved_and_candidates,sample_to_udnid){
  if(grepl("^e",myvar)){my_expspl="expression"}else{my_expspl="splicing"}
  
  if(my_expspl=="expression"){
    zs=get(myvar)%>%mutate(udnid=sample_to_udnid[sample_id])
  }else{
    zs=get(myvar)%>%mutate(udnid=sample_to_udnid[sample])
  }
  
  mytiss=sapply(strsplit(myvar,"_"),"[[",2)
  mybuild=sapply(strsplit(myvar,"_"),"[[",3)
  zs_solvedcases<-zs%>%dplyr::filter(udnid%in%solved_and_candidates$UDNID)
  solved_df<-data.frame()
  for(i in 1:nrow(solved_and_candidates)){
    this_udn_id=as.character(solved_and_candidates[i,"UDNID"]);this_ensg=as.character(solved_and_candidates[i,"ENSG"])
    found_gene_sample=zs_solvedcases%>%dplyr::filter(udnid==this_udn_id & ensg==this_ensg)
    if(nrow(found_gene_sample)>0){
      if(my_expspl=="expression"){
        solved_df=rbind(solved_df,
                        cbind(as.data.frame(found_gene_sample[,c("sample_id","ensg","zscore","rank_ind","status")]),solved_and_candidates[i,c("SOLVE_STATUS")],
                              tissue=mytiss,build=mybuild,exp_spl=my_expspl))
      }else{
        solved_df=rbind(solved_df,
                        cbind(as.data.frame(found_gene_sample[,c("sample","ensg","z","rank_ind","status")]),solved_and_candidates[i,c("SOLVE_STATUS")],
                              tissue=mytiss,build=mybuild,exp_spl=my_expspl))
      }
      
    }
  }
  if(my_expspl=="expression"){
    return(solved_df)
  }
  else{
    return(solved_df%>%rename(sample_id=sample,zscore=z))
  }
  
}
solved_and_candidates_zs<-data.frame()
for(this_zscore in myzscores){
  solved_and_candidates_zs=rbind(solved_and_candidates_zs,
                                 get_solved_samples(this_zscore,solved_and_candidates,sample_to_udnid))
}
# write_tsv(solved_and_candidates_zs,file="/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ZscoreDifferences/solved_outliers.tsv")
solved_and_candidates_zs=unique(read_tsv("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ZscoreDifferences/solved_outliers.tsv"))
solved_and_candidates_zs_wide=solved_and_candidates_zs%>%pivot_wider(names_from=build,values_from=zscore,id_cols = c(sample_id,ensg,status,SOLVE_STATUS,tissue,exp_spl))
# solved_hg19_hg38=solved_and_candidates_zs%>%dplyr::filter(build=="hg19" | build=="hg38")
# solved_hg38_chm13=solved_and_candidates_zs%>%dplyr::filter(build=="hg38" | build=="chm13ensembl")
ggplot(solved_and_candidates_zs_wide,aes(x=chm13ensembl,y=hg38,
                                         shape=exp_spl,color=chm13ensembl-hg38))+
  geom_point(aes(size=3,alpha=0.3))+
  # xlab("chm13")+
  
  stat_cor(label.y = 6,inherit.aes = F,aes(x=chm13ensembl,y=hg38),size=5)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 5,output.type="text",inherit.aes = F,aes(x=chm13ensembl,y=hg38),size=5) +
  scale_color_gradientn(colors=c("#dea700","#006B2F","#01427a"),limits = c(-1,1))+
  xlab("chm13")+
  # stat_cor(label.y = 6,inherit.aes = F,aes(x=hg19,y=hg38),size=5)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  # stat_regline_equation(label.y = 5,output.type="text",inherit.aes = F,aes(x=hg19,y=hg38),size=5) +
  # scale_color_gradientn(colors=c("#8C1515","#d07c10","#dea700"),limits = c(-1,1))+
  ggtitle("solved and candidate cases")+
  theme_bw(base_size = 20)

solved_and_candidates_zs_wide=solved_and_candidates_zs_wide%>%mutate(hg38hg19diff=hg38-hg19)%>%mutate(chm13hg38diff=chm13ensembl-hg38)%>%mutate(udn_id=sample_to_udnid[sample_id])

######OMIM
hg38_chm13_spl=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/zscores_genesdiff_hg38_chm13.soutliers.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))
colnames(hg38_chm13_spl)[c(4,7)]<-c("z_hg38","z_chm13")
hg19_hg38_spl=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/zscores_genesdiff_hg19_hg38.soutliers.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))
colnames(hg19_hg38_spl)[c(4,7)]<-c("z_hg19","z_hg38")

hg38_chm13_fibroblast_exp=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg38.chm13ensembl.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))%>%mutate(tissue="Fibroblast")%>%mutate(rank_build1=rank(((z_hg38)),ties.method="first"))%>%mutate(rank_build2=rank(((z_chm13)),ties.method="first"))%>%mutate(comparison="hg38:chm13")%>%dplyr::filter(abs(z_hg38)>.5 | abs(z_chm13)>.5)
hg38_chm13_blood_exp=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg38.chm13ensembl.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))%>%mutate(tissue="Blood")%>%mutate(rank_build1=rank(((z_hg38)),ties.method="first"))%>%mutate(rank_build2=rank(((z_chm13)),ties.method="first"))%>%mutate(comparison="hg38:chm13")%>%dplyr::filter(abs(z_hg38)>.5 | abs(z_chm13)>.5)
hg19_hg38_fibroblast_exp=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))%>%mutate(tissue="Fibroblast")%>%mutate(rank_build1=rank(((z_hg19)),ties.method="first"))%>%mutate(rank_build2=rank(((z_hg38)),ties.method="first"))%>%mutate(comparison="hg19:hg38")%>%dplyr::filter(abs(z_hg38)>.5 | abs(z_hg19)>.5)
hg19_hg38_blood_exp=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))%>%mutate(tissue="Blood")%>%mutate(rank_build1=rank(((z_hg19)),ties.method="first"))%>%mutate(rank_build2=rank(((z_hg38)),ties.method="first"))%>%mutate(comparison="hg19:hg38")%>%dplyr::filter(abs(z_hg38)>.5 | abs(z_hg19)>.5)


return_spl_exp<-function(comp,sample_to_tiss){
  comp=comp%>%dplyr::select(-transcript.x,-transcript.y,-chr.x)
  comp$tissue<-sample_to_tiss[comp$sample_id]
  comp_blood=comp%>%dplyr::filter(tissue=="Blood")
  comp_fibroblast=comp%>%dplyr::filter(tissue=="Fibroblast")
  return(list(comp_fibroblast,comp_blood))
}

hg38_chm13_spl=return_spl_exp(hg38_chm13_spl,sample_to_tiss)
hg38_chm13_fibroblast_spl=hg38_chm13_spl[[1]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg38:chm13")%>%dplyr::filter(z_chm13!=0 | z_hg38!=0)
hg38_chm13_blood_spl=hg38_chm13_spl[[2]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg38:chm13")%>%dplyr::filter(z_chm13!=0 | z_hg38!=0)
hg19_hg38_spl=return_spl_exp(hg19_hg38_spl,sample_to_tiss)
hg19_hg38_fibroblast_spl=hg19_hg38_spl[[1]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(comparison="hg19:hg38")%>%dplyr::filter(z_hg19!=0 | z_hg38!=0)
hg19_hg38_blood_spl=hg19_hg38_spl[[2]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(comparison="hg19:hg38")%>%dplyr::filter(z_hg19!=0 | z_hg38!=0)

solved_and_candidates<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_AND_candidate_dec2022.txt")
solved_and_candidates=solved_and_candidates%>%mutate(single_col=paste0(RDID,ENSG))

zs_hg19_hg38=rbind(hg19_hg38_blood_spl,hg19_hg38_fibroblast_spl,hg19_hg38_blood_exp,hg19_hg38_fibroblast_exp)%>%dplyr::mutate(is_omim=ifelse(gene%in%omim_genes,"omim","non-omim"),status=sample_to_status[sample_id]) %>%mutate(single_col=paste0(sample_id,gene)) #%>%select(-z_hg19,-z_hg38)
zs_hg38_chm13=rbind(hg38_chm13_blood_spl,hg38_chm13_fibroblast_spl,hg38_chm13_fibroblast_exp,hg38_chm13_blood_exp)%>%dplyr::mutate(is_omim=ifelse(gene%in%omim_genes,"omim","non-omim"),status=sample_to_status[sample_id])%>%mutate(single_col=paste0(sample_id,gene)) #%>%select(-z_chm13,-z_hg38))

zs_hg19_hg38_solved=zs_hg19_hg38%>%dplyr::filter(single_col%in%solved_and_candidates$single_col)
zs_hg38_chm13_solved=zs_hg38_chm13%>%dplyr::filter(single_col%in%solved_and_candidates$single_col)

zs_hg19_hg38_diffs=zs_hg19_hg38%>%dplyr::filter((abs(z_hg19)>3 | abs(z_hg38)>3) & abs(z_diff)>2 & status=="Case"&is_omim=="omim")%>%dplyr::filter(rank_build1<30 | rank_build2<30)
zs_hg38_chm13_diffs=zs_hg38_chm13%>%dplyr::filter((abs(z_chm13)>3 | abs(z_hg38)>3) & abs(z_diff)>2&status=="Case"&is_omim=="omim") %>%dplyr::filter(rank_build1<30 | rank_build2<30)


#######pull the OMIM genes w large zdiff (z>3, z<1)
get_omim_list<-function(z,outlier_threshold=3,build_reference,build_to_compare,tissue,omim_genes,chrY_chm13_genes){
  if(build_reference=="z_chm13" | build_to_compare=="z_chm13"){
    print("chm13")
    print(nrow(z))
    z=z%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))
    print(nrow(z))
  }
  build_reference_name=gsub("z_","",build_reference)
  build_compare_name=gsub("z_","",build_to_compare)
  filt_zs=z%>%dplyr::filter(get(build_reference)>=outlier_threshold & get(build_to_compare)<1)
  in_omim=filt_zs%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))%>%dplyr::mutate(is_omim=ifelse(gene%in%omim_genes,"omim","non-omim")) %>%
    mutate(rank_diff=abs(rank_build1-rank_build2))%>%
    group_by(gene,comparison, tissue)%>%summarise(max_zdiff=max(abs(z_diff)),max_rankdiff=max(rank_diff))%>%
    mutate(comparison_direction=paste0(build_reference_name,":",build_compare_name))
  #%>%pull(gene)%>%unique()%>%length()
  if(nrow(in_omim)>0){
    return(as.data.frame(in_omim))
  }
  
}
eOutlier_omim_diff_df=data.frame()
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg38_chm13_blood_exp,3,"z_hg38","z_chm13","blood",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg38_chm13_fibroblast_exp,3,"z_hg38","z_chm13","fibroblast",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg19_hg38_blood_exp,3,"z_hg19","z_hg38","blood",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg19_hg38_fibroblast_exp,3,"z_hg19","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg38_chm13_blood_exp,3,"z_chm13","z_hg38","blood",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg38_chm13_fibroblast_exp,3,"z_chm13","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg19_hg38_blood_exp,3,"z_hg38","z_hg19","blood",omim_genes,chrY_chm13_genes))
eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list(hg19_hg38_fibroblast_exp,3,"z_hg38","z_hg19","fibroblast",omim_genes,chrY_chm13_genes))

#overlap w solved
eOutlier_omim_diff_solvedgene=eOutlier_omim_diff_df%>%dplyr::filter(gene%in%solved_and_candidates$ENSG)
hg38_chm13_blood_spl%>%dplyr::filter(gene%in%omim_diff_solvedgene$gene)%>%dplyr::filter(abs(z_diff)>2)%>%mutate(status=sample_to_status[sample_id])

sOutlier_omim_diff_df=data.frame()
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg38_chm13_blood_spl,3,"z_hg38","z_chm13","blood",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg38_chm13_fibroblast_spl,3,"z_hg38","z_chm13","fibroblast",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg19_hg38_blood_spl,3,"z_hg19","z_hg38","blood",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg19_hg38_fibroblast_spl,3,"z_hg19","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg38_chm13_blood_spl,3,"z_chm13","z_hg38","blood",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg38_chm13_fibroblast_spl,3,"z_chm13","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg19_hg38_blood_spl,3,"z_hg38","z_hg19","blood",omim_genes,chrY_chm13_genes))
sOutlier_omim_diff_df=rbind(sOutlier_omim_diff_df,get_omim_list(hg19_hg38_fibroblast_spl,3,"z_hg38","z_hg19","fibroblast",omim_genes,chrY_chm13_genes))

#overlap w solved
solved_and_candidates<-fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_AND_candidate_dec2022.txt")
omim_diff_solvedgene=sOutlier_omim_diff_df%>%dplyr::filter(gene%in%solved_and_candidates$ENSG)
hg38_chm13_blood_spl%>%dplyr::filter(gene%in%omim_diff_solvedgene$gene)%>%dplyr::filter(abs(z_diff)>2)%>%mutate(status=sample_to_status[sample_id])

#### odds ratios
get_OR<-function(fibroblast_hg38eOutlier,this_splexp,this_tissue,sample_to_status=sample_to_status){
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

OR_fibroblast_exp_chm13=get_OR(eoutliers_fibroblast_chm13,"expression","fibroblast")
OR_blood_exp_chm13=get_OR(eoutliers_blood_chm13,"expression","blood")
OR_fibroblast_spl_chm13=get_OR(soutliers_fibroblast_chm13,"splicing","fibroblast")
OR_blood_spl_chm13=get_OR(soutliers_blood_chm13,"splicing","blood")
OR_fibroblast_exp_hg38=get_OR(eoutliers_fibroblast_hg38,"expression","fibroblast")
OR_blood_exp_hg38=get_OR(eoutliers_blood_hg38,"expression","blood")
OR_fibroblast_spl_hg38=get_OR(soutliers_fibroblast_hg38,"splicing","fibroblast")
OR_blood_spl_hg38=get_OR(soutliers_blood_hg38,"splicing","blood")
OR_fibroblast_exp_hg19=get_OR(eoutliers_fibroblast_hg19,"expression","fibroblast")
OR_blood_exp_hg19=get_OR(eoutliers_blood_hg19,"expression","blood")
OR_fibroblast_spl_hg19=get_OR(soutliers_fibroblast_hg19,"splicing","fibroblast")
OR_blood_spl_hg19=get_OR(soutliers_blood_hg19,"splicing","blood")
all_OR=rbind(cbind(rbind(OR_fibroblast_exp_chm13,OR_blood_exp_chm13,OR_fibroblast_spl_chm13,OR_blood_spl_chm13),build="chm13",uniq="not unique"),
             cbind(rbind(OR_fibroblast_exp_hg38,OR_blood_exp_hg38,OR_fibroblast_spl_hg38,OR_blood_spl_hg38),build="hg38",uniq="not unique"),
             cbind(rbind(OR_fibroblast_exp_hg19,OR_blood_exp_hg19,OR_fibroblast_spl_hg19,OR_blood_spl_hg19),build="hg19",uniq="not unique"))
ggplot(all_OR%>%dplyr::filter(tissue=="blood"),
       aes(x=spl_or_exp,y=`odds_ratio.odds ratio`,color=build))+
  geom_point(aes(size = 2),position=position_dodge(width=.8)) +
  scale_color_manual(values=c("#8C1515","#dea700","#01427a"))+
  geom_hline(yintercept=1,color="#888888")+
  theme_bw()+
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),width=0.4, position = position_dodge(width=.8))+
  # labs(title = 'Odds an outlier is in a case as compared to a control', x = '', y = 'OR')
  ggtitle("Odds an outlier is in a case as compared to a control")+xlab("")+ylab("OR")
