library("biomaRt")
library(data.table)
library(dplyr)
library(ggpubr)


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
  hgnc_more=as.data.frame(fread("/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/hgnc_convert.txt"))
  ensg_to_hgnc_more=hgnc_more[,10]
  names(ensg_to_hgnc_more)=(hgnc_more[,2])
  combined_ensg_to_hgnc=unique(rbind(stack(ensg_to_hgnc_more),stack(hgnc_to_ensg)))
  colnames(combined_ensg_to_hgnc)<-c("ensg","hgnc")
  # write_tsv(combined_ensg_to_hgnc, "/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/hgnc_ensg_convert_complete.txt")
  ref_dir="/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/"
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

md=fread("/oak/stanford/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv")
sample_to_status=md$affected_status
names(sample_to_status)=md$sample_id
sample_to_tiss=md$source_of_RNA
names(sample_to_tiss)=md$sample_id

chrY_chm13_genes=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/chm13_chrY.txt",header = F)
colnames(chrY_chm13_genes)<-c("hgnc","ensembl","chm13_id")
chrY_chm13_genes<-chrY_chm13_genes%>%mutate(gene=gsub("\\..*","",ensembl))


#########SUMMARY
###SUMMARY STATS
sOutliers_hg19<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg19.sorted_dedupOptical_minMQ255.bed.gz")
sOutliers_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg38.sorted_dedupOptical_minMQ255.bed.gz")
sOutliers_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.chm13.sorted_dedupOptical_minMQ255.bed.gz")
sOutliers_chm13ensembl<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/Splicing/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.chm13.sorted_dedupOptical_minMQ255.bed.gz")
colnames(sOutliers_hg19)<-c("chr","pos1","pos2","clus","sample","z","enst","gene","x")
colnames(sOutliers_hg38)<-c("chr","pos1","pos2","clus","sample","z","junc","junc_support","strand","splicesite","a","b","c","noveltype","novel_donor","novel_acceptor","enst","gene","x")
colnames(sOutliers_chm13)<-c("chr","pos1","pos2","clus","sample","z","enst","gene","chm13_id","gene2","x")
colnames(sOutliers_chm13ensembl)<-c("chr","pos1","pos2","clus","sample","z","enst","gene","chm13_id","gene2","x")
sOutliers_hg19=sOutliers_hg19%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])
sOutliers_chm13=sOutliers_chm13%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])
sOutliers_hg38=sOutliers_hg38%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])
sOutliers_chm13ensembl=sOutliers_chm13ensembl%>%mutate(tissue=sample_to_tiss[sample])%>%mutate(status=sample_to_status[sample])
sOutliers_chm13=sOutliers_chm13%>%group_by(sample)%>%mutate(ind_rank=rank(desc(abs(z)),ties.method="first"))
sOutliers_hg38=sOutliers_hg38%>%group_by(sample)%>%mutate(ind_rank=rank(desc(abs(z)),ties.method="first"))
sOutliers_hg19=sOutliers_hg19%>%group_by(sample)%>%mutate(ind_rank=rank(desc(abs(z)),ties.method="first"))
sOutliers_hg19_fibroblast=sOutliers_hg19%>%dplyr::filter(tissue=="Fibroblast")
sOutliers_hg19_blood=sOutliers_hg19%>%dplyr::filter(tissue=="Blood")
sOutliers_hg38_fibroblast=sOutliers_hg38%>%dplyr::filter(tissue=="Fibroblast")
sOutliers_hg38_blood=sOutliers_hg38%>%dplyr::filter(tissue=="Blood")
sOutliers_chm13ensembl_fibroblast=sOutliers_chm13ensembl%>%dplyr::filter(tissue=="Fibroblast")
sOutliers_chm13ensembl_blood=sOutliers_chm13ensembl%>%dplyr::filter(tissue=="Blood")
sOutliers_chm13_fibroblast=sOutliers_chm13%>%dplyr::filter(tissue=="Fibroblast")
sOutliers_chm13_blood=sOutliers_chm13%>%dplyr::filter(tissue=="Blood")


myzscores=c("sOutliers_hg19_fibroblast","sOutliers_hg19_blood",
            "sOutliers_hg38_fibroblast","sOutliers_hg38_blood",
            "sOutliers_chm13_fibroblast","sOutliers_chm13_blood")
summ_stats_spl_df=data.frame()
all_soutliers=data.frame()
for(this_zscore in myzscores){
  my_build=sapply(str_split(this_zscore,"_"),"[[",2)
  outliers=get(this_zscore)%>%dplyr::filter(abs(z)>3)
  all_soutliers=rbind(all_soutliers,cbind(outliers, build=my_build,exp_spl="sOutlier"))
  num_outliers=outliers%>%nrow()
  num_outlier_genes=length(unique(outliers$gene))
  uniq_omim_gene_overlap=outliers%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene %in% omim_genes)%>%pull(gene)%>%unique()%>%length()
  summ_stats_spl_df=rbind(as.data.frame(summ_stats_spl_df),
                          as.data.frame(cbind(this_zscore,num_outliers,num_outlier_genes,uniq_omim_gene_overlap, exp_spl="sOutlier",build=my_build)))
}

num_outliers_per_sample=rbind(all_soutliers%>%dplyr::filter(z>3)%>%group_by(sample,exp_spl,status,build,tissue)%>%summarise(num_outliers=n())%>%rename(sample_id=sample),
                              all_eoutliers%>%dplyr::filter(zscore>3)%>%group_by(sample_id,exp_spl,status,build,tissue)%>%summarise(num_outliers=n()),
                              all_eoutliers%>%dplyr::filter(zscore<(-3))%>%group_by(sample_id,exp_spl,status,build,tissue)%>%summarise(num_outliers=n()))
# num_outliers_per_sample=num_outliers_per_sample[,c(1:3,5,7)]
# colnames(num_outliers_per_sample)[3:5]<-c("sOutlier","Under eOutlier","Over eOutlier")
# num_outliers_per_sample_melted=melt(num_outliers_per_sample)


summstats_spl=summ_stats_spl_df%>%separate(this_zscore, c("splinfo","build","tissue"),sep = "_")%>%mutate(exp_spl="spl")%>%dplyr::select(-splinfo)
summstats_exp=summ_stats_exp_df%>%separate(this_zscore, c("tissue","build"),sep = "_")%>%mutate(exp_spl="exp")
summstats_all=rbind(summstats_exp,summstats_spl) %>%dplyr::filter(tissue=="blood")
ggplot(num_outliers_per_sample,aes(x=exp_spl,y=as.numeric(as.character(uniq_omim_gene_overlap)),
                         fill=build))+
  geom_bar(position="dodge",stat="identity")+
  facet_wrap(~tissue)+
  xlab("")+ylab("Number of outliers")+
  geom_text(aes(label=num_outliers), position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()
melted_summstats_all=melt(setDT(summstats_all), id.vars = c("tissue","build","exp_spl"), variable.name = "summ_stat")
melted_summstats_all$value<-as.numeric(melted_summstats_all$value)
#fig 3a
ggplot(melted_summstats_all,aes(x=exp_spl,y=value,
                         group=build,fill=summ_stat))+
  geom_bar(position="dodge",stat="identity")+
  facet_wrap(~tissue,scales = "free_y")+
  xlab("")+ylab("Number of OMIM outliers")+
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()

######COMPARISONS
hg38_chm13=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/zscores_genesdiff_hg38_chm13.soutliers.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))
colnames(hg38_chm13)[c(4,7)]<-c("z_hg38","z_chm13")
hg19_hg38=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/zscores_genesdiff_hg19_hg38.soutliers.txt.gz")%>%mutate(gene=gsub("\\..*","",gene))
colnames(hg19_hg38)[c(4,7)]<-c("z_hg19","z_hg38")

return_spl_exp<-function(comp,sample_to_tiss){
  comp=comp%>%dplyr::select(-transcript.x,-transcript.y,-chr.x)
  comp$tissue<-sample_to_tiss[comp$sample_id]
  comp_blood=comp%>%dplyr::filter(tissue=="Blood")
  comp_fibroblast=comp%>%dplyr::filter(tissue=="Fibroblast")
  return(list(comp_fibroblast,comp_blood))
}

hg38_chm13_spl=return_spl_exp(hg38_chm13,sample_to_tiss)
hg38_chm13_fibroblast_spl=hg38_chm13_spl[[1]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg38:chm13")
hg38_chm13_blood_spl=hg38_chm13_spl[[2]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg38:chm13")
hg19_hg38_spl=return_spl_exp(hg19_hg38,sample_to_tiss)
hg19_hg38_fibroblast_spl=hg19_hg38_spl[[1]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(comparison="hg19:hg38")
hg19_hg38_blood_spl=hg19_hg38_spl[[2]] %>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%
  mutate(rank_build2=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(comparison="hg19:hg38")

get_outlier_prop<-function(z,outlier_threshold=3,build_reference,build_to_compare,tissue,omim_genes,chrY_chm13_genes){
  if(build_reference=="z_chm13" | build_to_compare=="z_chm13"){
    print("chm13")
    print(nrow(z))
    z=z%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))
    print(nrow(z))
  }
  print(paste0("What proportion of ",build_reference, " outliers are consistent in ",build_to_compare,"?"))
  num_outliers_consistent=round(z%>%dplyr::filter(get(build_reference)>=outlier_threshold & get(build_to_compare)>=outlier_threshold)%>%nrow()/
                                  z%>%dplyr::filter(get(build_reference)>=outlier_threshold)%>%nrow()*100,outlier_threshold)
  print(paste0("percent of outliers that are consistent: ",num_outliers_consistent,"%"))
  num_outliers_notseenother=round(z%>%dplyr::filter(get(build_reference)>=outlier_threshold & is.na(get(build_to_compare)))%>%nrow()/
                                    z%>%dplyr::filter(get(build_reference)>=outlier_threshold)%>%nrow()*100,outlier_threshold)
  print(paste0("percent not seen in the other build: ",num_outliers_notseenother ,"%"))
  # for_df<-c(paste0("prop of ",build_reference," (reference) outliers are also outliers in ",build_to_compare),tissue,
  for_df<-c(paste0(build_reference,"-",build_to_compare),tissue,
            num_outliers_consistent,
            num_outliers_notseenother)
  for(i in c(1,2,3)){
    denom=z%>%dplyr::filter(get(build_reference)>=outlier_threshold)%>%nrow()
    filt_zs=z%>%dplyr::filter(get(build_reference)>=outlier_threshold & get(build_to_compare)<i)
    # in_omim=filt_zs%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))%>%dplyr::filter(gene%in%omim_genes)%>%pull(gene)%>%unique()%>%length()
    numerator=filt_zs%>%nrow()
    frac=numerator/denom*100
    print(paste0("z<",i, " is ", round(frac,outlier_threshold), "%"))
    for_df<-c(for_df,round(frac,3)) #,in_omim)
    # for_df<-c(for_df,round(frac,3),in_omim)
  }
  return_df=as.data.frame(t(for_df))
  # colnames(return_df)<-c("comparison","tissue","consistent","not_in_other_build","z<1", "omim_z<1",
  #                        "z<2", "omim_z<2","z<3", "omim_z<3")
  colnames(return_df)<-c("comparison","tissue","consistent","not_in_other_build","z<1",
                         "z<2","z<3")
  return(return_df)
}
get_soutlier_consistency_df<-function(){

  sOutlier_consistency_df=data.frame()
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg38_chm13_blood_spl,3,"z_hg38","z_chm13","blood",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg38_chm13_fibroblast_spl,3,"z_hg38","z_chm13","fibroblast",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg19_hg38_blood_spl,3,"z_hg19","z_hg38","blood",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg19_hg38_fibroblast_spl,3,"z_hg19","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg38_chm13_blood_spl,3,"z_chm13","z_hg38","blood",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg38_chm13_fibroblast_spl,3,"z_chm13","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg19_hg38_blood_spl,3,"z_hg38","z_hg19","blood",omim_genes,chrY_chm13_genes))
  sOutlier_consistency_df=rbind(sOutlier_consistency_df,get_outlier_prop(hg19_hg38_fibroblast_spl,3,"z_hg38","z_hg19","fibroblast",omim_genes,chrY_chm13_genes))
  return(sOutlier_consistency_df)
}
sOutlier_consistency_df=get_soutlier_consistency_df()


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
    in_omim=filt_zs%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))%>%dplyr::filter(gene%in%omim_genes) %>%
      mutate(rank_diff=abs(rank_build1-rank_build2))%>%
      group_by(gene,comparison, tissue)%>%summarise(max_zdiff=max(abs(z_diff)),max_rankdiff=max(rank_diff))%>%
      mutate(comparison_direction=paste0(build_reference_name,":",build_compare_name))
      #%>%pull(gene)%>%unique()%>%length()
    if(nrow(in_omim)>0){
      return(as.data.frame(in_omim))
    }
  
}
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
solved_and_candidates<-fread("/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_AND_candidate_dec2022.txt")
omim_diff_solvedgene=sOutlier_omim_diff_df%>%dplyr::filter(gene%in%solved_and_candidates$ENSG)
hg38_chm13_blood_spl%>%dplyr::filter(gene%in%omim_diff_solvedgene$gene)%>%dplyr::filter(abs(z_diff)>2)%>%mutate(status=sample_to_status[sample_id])

####rank comparison
spl_all_comparisons=rbind(hg19_hg38_blood_spl,hg19_hg38_fibroblast_spl,
                      hg38_chm13_blood_spl,hg38_chm13_fibroblast_spl)
spl_all_comparisons_top20=spl_all_comparisons%>%dplyr::filter(rank_build1<20 | rank_build2<20)

ggplot(spl_all_comparisons_top20,
       aes(x=rank_build1,y=rank_build2,alpha=0.3,color=tissue))+
  geom_point()+
  scale_x_continuous(trans="log10")+
  theme_bw()+
  facet_wrap(~comparison)+
  stat_regline_equation(label.x=3,aes(label = ..rr.label..)) +
  ggtitle("Splicing Outlier Rank Comparison")

## diff exp overlap
diffexp_blood_hg38_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/Blood.dedupOptical_minMQ255.diff_expression_chm13ensembl_vs_hg38.limma_voom_dream.significant.txt.gz")%>%dplyr::filter(abs(logFC)>1)
diffexp_blood_hg19_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/Blood.dedupOptical_minMQ255.diff_expression_hg38_vs_hg19.limma_voom_dream.significant.txt.gz")%>%dplyr::filter(abs(logFC)>1)
diffexp_blood_hg19_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/Blood.dedupOptical_minMQ255.diff_expression_chm13ensembl_vs_hg19.limma_voom_dream.significant.txt.gz")%>%dplyr::filter(abs(logFC)>1)

diffexp_fibroblast_hg38_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/Fibroblast.dedupOptical_minMQ255.diff_expression_chm13ensembl_vs_hg38.limma_voom_dream.significant.txt.gz")%>%dplyr::filter(abs(logFC)>1)
diffexp_fibroblast_hg19_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/Fibroblast.dedupOptical_minMQ255.diff_expression_hg38_vs_hg19.limma_voom_dream.significant.txt.gz")%>%dplyr::filter(abs(logFC)>1)
diffexp_fibroblast_hg19_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/Fibroblast.dedupOptical_minMQ255.diff_expression_chm13ensembl_vs_hg19.limma_voom_dream.significant.txt.gz")%>%dplyr::filter(abs(logFC)>1)
compare_outliers_diffspl<-function(z,comp1,comp2,diff_exp,omim_genes,tissue,comparison){
  diff_exp=diff_exp%>%dplyr::filter(abs(logFC)>1 & adj.P.Val<0.05)
  z_diff_exp<-z%>%mutate(diff_exp=ifelse(gene %in% diff_exp$gene,"differentially_quantified","not_differentially_quantified"))
  z_diff_exp_outliers<-z_diff_exp%>%dplyr::filter(abs(get(comp1))>3 | abs(get(comp2))>3)
  z_diff_exp_largediff<-z_diff_exp%>%dplyr::filter((abs(get(comp1))>3 & abs(get(comp2))<1 ) |
                                                     (abs(get(comp1))<1 & abs(get(comp2))>3 )
  )
  uniq_genes_quantorno=z_diff_exp_largediff%>%ungroup()%>%dplyr::select(gene,diff_exp)%>%unique()
  mytable=table(uniq_genes_quantorno$diff_exp)
  mydf=cbind(comparison=comparison,tissue=tissue,
             total_diffexp_genes=nrow(diff_exp),
             diffexp_genes_no_outliers=table(diff_exp$gene %in% z_diff_exp_outliers$gene)["FALSE"],
             diffexp_gene_yes_outliers=table(diff_exp$gene %in% z_diff_exp_outliers$gene)["TRUE"],
             uniq_gene_outliers=length(unique(z_diff_exp_outliers$gene)),
             num_large_differences=length(unique(z_diff_exp_largediff$gene)),
             t(data.frame(unclass(mytable))))
  
  # uniq_omim_gene_overlap=outliers%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene %in% omim_genes)%>%pull(gene)%>%unique()%>%length()
  return(as.data.frame(mydf))
}

outlier_diffspl_table=rbind(compare_outliers_diffspl(hg19_hg38_blood_spl,"z_hg19","z_hg38",diffexp_blood_hg19_hg38,omim_genes,"blood","hg19:chm38"),
                            compare_outliers_diffspl(hg19_hg38_fibroblast_spl,"z_hg19","z_hg38",diffexp_fibroblast_hg19_hg38,omim_genes,"fibroblast","hg19:chm38"),
                            compare_outliers_diffspl(hg38_chm13_blood_spl,"z_hg38","z_chm13",diffexp_blood_hg38_chm13,omim_genes,"blood","hg38:chm13"),
                            compare_outliers_diffspl(hg38_chm13_fibroblast_spl,"z_hg38","z_chm13",diffexp_fibroblast_hg38_chm13,omim_genes,"fibroblast","hg38:chm13"))


## annotation specific
get_soutlier_ONEannoONLY<-function(soutliers,gene_id_filt,build,main_table,max_rank=20){
  soutliers=soutliers%>%dplyr::filter(ind_rank<=max_rank&status=="Case")

  if(build=="chm13"){
    sOutliers_annoONLY=soutliers%>%dplyr::filter(chm13_id%in%gene_id_filt)
    sOutliers_annoONLY_outliers=sOutliers_annoONLY%>%dplyr::filter(abs(z)>3)
    sOutliers_annoONLY_outliers_case=sOutliers_annoONLY_outliers%>%dplyr::filter(status=="Case")
    sOutliers_annoONLY_outliers_control=sOutliers_annoONLY_outliers%>%dplyr::filter(status=="Control")
    sOutliers_top20=sOutliers_annoONLY%>%dplyr::filter(ind_rank<=20&z>0)
    sOutliers_top20_caseonly=sOutliers_top20%>%dplyr::filter(status=="Case")
    print(paste0("num outliers generally: ",length(unique(sOutliers_annoONLY_outliers$chm13_id))))
    print(paste0("uniq top 20 case and z>0: ", length(unique(sOutliers_top20_caseonly$chm13_id))))
    
    print(paste0("num outliers cases: ",length(unique(sOutliers_annoONLY_outliers_case$chm13_id))))
    
    genetable_outliers_cases=main_table%>%dplyr::filter(gene_id%in%sOutliers_annoONLY_outliers_case$chm13_id)
    print(genetable_outliers_cases%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
            dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table())
    
    print(paste0("num outliers controls: ",length(unique(sOutliers_annoONLY_outliers_control$chm13_id))))
    genetable_outliers_controls=main_table%>%dplyr::filter(gene_id%in%sOutliers_annoONLY_outliers_control$chm13_id)
    print(genetable_outliers_controls%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
            dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table())
    
    
    nonissue_genes=main_table%>%dplyr::filter(gene_id%in%sOutliers_top20_caseonly$chm13_id)%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
      dplyr::filter(any_issue=='non_issue')
    print(table(nonissue_genes$any_clin_relevant))
    
    print('clinically relevant?')
    print((main_table%>%dplyr::filter(gene_id %in% sOutliers_annoONLY_outliers$gene & any_clin_relevant=='clinically_relevant')%>%dplyr::select(gene_id,gene_name)))
    print(table(main_table%>%dplyr::filter(gene_id %in% sOutliers_annoONLY_outliers$gene )%>%dplyr::select(tissue,any_clin_relevant)))
    
    return(main_table%>%dplyr::filter(gene_id%in%sOutliers_top20_caseonly$chm13_id))
    
    
  }else{
    sOutliers_annoONLY=soutliers%>%dplyr::filter(gene%in%gene_id_filt)
    sOutliers_annoONLY_outliers=sOutliers_annoONLY%>%dplyr::filter(abs(z)>3)
    sOutliers_annoONLY_outliers_case=sOutliers_annoONLY_outliers%>%dplyr::filter(status=="Case")
    sOutliers_annoONLY_outliers_control=sOutliers_annoONLY_outliers%>%dplyr::filter(status=="Control")
    sOutliers_top20=sOutliers_annoONLY%>%dplyr::filter(ind_rank<=20&z>3)
    sOutliers_top20_caseonly=sOutliers_top20%>%dplyr::filter(status=="Case")
    print(paste0("num outliers generally: ",length(unique(sOutliers_annoONLY_outliers$gene))))
    print(paste0("uniq top 20 case and z>3: ", length(unique(sOutliers_top20_caseonly$gene))))
  
    
    print(paste0("num outliers cases: ",length(unique(sOutliers_annoONLY_outliers_case$gene))))
    genetable_outliers_cases=main_table%>%dplyr::filter(gene_id%in%sOutliers_annoONLY_outliers_case$gene)
    print(genetable_outliers_cases%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
            dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table())
    
    print(paste0("num outliers controls: ",length(unique(sOutliers_annoONLY_outliers_control$gene))))
    genetable_outliers_controls=main_table%>%dplyr::filter(gene_id%in%sOutliers_annoONLY_outliers_control$gene)
    print(genetable_outliers_controls%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
            dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table())
    
    
    nonissue_genes=main_table%>%dplyr::filter(gene_id%in%sOutliers_top20_caseonly$gene)%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
      dplyr::filter(any_issue=='non_issue')
    print(table(nonissue_genes$any_clin_relevant))
    print(nonissue_genes%>%dplyr::filter(!is.na(any_clin_relevant)))
    
    print((main_table%>%dplyr::filter(gene_id %in% sOutliers_annoONLY_outliers$gene & any_clin_relevant=='clinically_relevant')%>%dplyr::select(gene_id,gene_name)))
    print(table(main_table%>%dplyr::filter(gene_id %in% sOutliers_annoONLY_outliers$gene )%>%dplyr::select(tissue,any_clin_relevant)))
    
    return(main_table%>%dplyr::filter(gene_id%in%sOutliers_top20_caseonly$gene))
  }

  
}
chm13_only=main_table%>%dplyr::filter(hg38_chm13_annotation_comparison=='chm13_specific'|hg38_chm13_annotation_comparison=='chm13_spec_paralog')
hg38vschm13_only=main_table%>%dplyr::filter(hg38_chm13_annotation_comparison=='hg38_specific'|hg38_chm13_annotation_comparison=='hg38_spec_inferred'|hg38_chm13_annotation_comparison=='hg38_spec_gene_collapsed_in_chm13')
hg38vshg19_only=main_table%>%dplyr::filter(hg19_hg38_annotation_comparison=='hg38_specific')
hg19_only=main_table%>%dplyr::filter(hg19_hg38_annotation_comparison=='hg19_specific')
sOutliers_hg19=sOutliers_hg19%>%mutate(gene=gsub("\\..*","",gene))
sOutliers_hg38=sOutliers_hg38%>%mutate(gene=gsub("\\..*","",gene))
sOutliers_chm13=sOutliers_chm13%>%mutate(gene=gsub("\\..*","",gene))

hg19_spec_soutliers=sOutliers_hg19%>%dplyr::filter(gene%in%hg19_spec_genes&ind_rank<=20&status=="Case")%>%mutate(unique_in="hg19")%>%mutate(comparison="hg19:hg38")
hg38compared2hg19_spec_soutliers=sOutliers_hg38%>%dplyr::filter(gene%in%hg38_compared2hg19_spec_genes&ind_rank<=20&status=="Case")%>%mutate(unique_in="hg38_ascomparedtohg19")%>%dplyr::filter(build=='hg38')%>%mutate(comparison="hg19:hg38")
hg38compared2chm13_spec_soutliers=sOutliers_hg38%>%dplyr::filter(gene%in%hg38_compared2chm13_spec_genes&ind_rank<=20&status=="Case")%>%mutate(unique_in="hg38_ascomparedtochm13")%>%dplyr::filter(build=='hg38')%>%mutate(comparison="hg38:chm13")
chm13_spec_soutliers=sOutliers_chm13%>%dplyr::filter(gene%in%chm13_spec_genes&ind_rank<=20&status=="Case")%>%mutate(unique_in="chm13")%>%mutate(comparison="hg38:chm13")

all_spec_soutliers_top20=rbind(hg19_spec_soutliers,hg38compared2hg19_spec_soutliers,hg38compared2chm13_spec_soutliers,chm13_spec_soutliers)

hg19_top20_maintable=get_soutlier_ONEannoONLY(sOutliers_hg19,hg19_only$gene_id,'hg19',main_table)
hg38vs19_top20_maintable=get_soutlier_ONEannoONLY(sOutliers_hg38,hg38vshg19_only$gene_id,'hg38',main_table)
hg38vschm13_top20_maintable=get_soutlier_ONEannoONLY(sOutliers_hg38,hg38vschm13_only$gene_id,'hg38',main_table)
chm13_top20_maintable=get_soutlier_ONEannoONLY(sOutliers_chm13,chm13_only$gene_id,'chm13',main_table)

hg19_top20_maintable%>%group_by(gene_id)%>%mutate(any_issue=ifelse(!is.na(hg19_issue)|!is.na(hg19_blacklist_reason),"issue","non_issue"))%>%dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table()
hg38vs19_top20_maintable%>%group_by(gene_id)%>%mutate(any_issue=ifelse(!is.na(hg38_issue)|!is.na(hg38_blacklist_reason),"issue","non_issue"))%>%dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table()
hg38vschm13_top20_maintable%>%group_by(gene_id)%>%mutate(any_issue=ifelse(!is.na(hg38_issue)|!is.na(hg38_blacklist_reason),"issue","non_issue"))%>%dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table()
chm13_top20_maintable%>%group_by(gene_id)%>%mutate(any_issue=ifelse(!is.na(chm13_issue)|!is.na(chm13_blacklist_reason),"issue","non_issue"))%>%dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table()
