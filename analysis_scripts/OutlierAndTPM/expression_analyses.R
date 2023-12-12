library("biomaRt")
library(data.table)
library(tidyverse)

chrY_chm13_genes=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/chm13_chrY.txt",header = F)
colnames(chrY_chm13_genes)<-c("hgnc","ensembl","chm13_id")
chrY_chm13_genes<-chrY_chm13_genes%>%mutate(gene=gsub("\\..*","",ensembl))

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
                           ifelse(confidence == 2, "linkage_mapping", 
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

opentargets=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/opentargetscores.txt")
colnames(opentargets)<-c("target_id","ensg","target_score","x")
opentargets_summ=opentargets%>%group_by(ensg)%>% dplyr::filter(target_score ==max(target_score)) %>%dplyr::select(-x) #

cancer_genes=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/cosmic_cancer_genes.tsv")
cancer_genes$ensg=gsub("\\..*","",str_extract(cancer_genes$Synonyms,"ENSG(.*)\\."))


md=fread("/oak/stanford/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv")
sample_to_status=md$affected_status
names(sample_to_status)=md$sample_id
sample_to_tiss=md$source_of_RNA
names(sample_to_tiss)=md$sample_id
#expression generally
filter_tpms = function(df, people_filt = .2, tpm_filt = 0.1,min_sample_size=30) {
  genes = unique(df$gene)
  ind_filt_both = round(people_filt*length(unique(df$sample))) #this is 20% of people
  exp_df=(df)%>%group_by(gene)%>%mutate(n = sum(tpm >=tpm_filt))%>%dplyr::filter(n>ind_filt_both)
  return(exp_df%>%dplyr::select(-n))
}
read_in_all_tpms<-function(tissues,builds,omim_genes,opentargets_summ){
  build_to_dir_dic=c("hg19","hg38","chm13")
  names(build_to_dir_dic)=c("hg19","hg38","chm13ensembl")
  full_df<-data.frame()
  for(this_tissue in tissues){
    for(this_build in builds){
      
      this_filename=paste0("/oak/stanford/groups/smontgom/shared/UDN/Output/",build_to_dir_dic[this_build],"/eOutliers/",this_tissue,".",this_build,".dedupOptical_minMQ255_rsem.genes.results")
      this_df=fread(this_filename)
      colnames(this_df)=c("sample","gene","transcript","tpm")
      percent_people_20=length(unique(this_df$sample))*.2
      threshold_num_people=ifelse(percent_people_20<30,30,percent_people_20)
      print(length(unique(this_df$gene)))
      this_df=this_df%>%group_by(gene)%>%mutate(to_filter=ifelse(sum(tpm>0.15)<threshold_num_people,"failed","pass"))%>%
        dplyr::filter(to_filter=="pass")%>%dplyr::select(-to_filter)%>%
        mutate(tissue=this_tissue,build=this_build)
      print(length(unique(this_df$gene)))
      full_df=rbind(full_df,this_df)
    }
  }
  full_df=full_df%>%mutate(gene_fullname=gene,gene=gsub("\\..*","",gene_fullname),omim=ifelse(gene %in% omim_genes,"omim","not_omim"))
  return(full_df)
}
all_tpms=read_in_all_tpms(tissues=c("Blood","Fibroblast","iPS","iPSC_NPC","Muscle","PBMC"),builds=c("hg19","hg38","chm13ensembl"),omim_genes,opentargets_summ)
exp_with_omim=all_tpms%>%dplyr::filter(omim=="omim")%>%pull(gene)%>%unique()
table(unique(all_tpms$gene)%in%opentargets_summ$ensg); length(unique(opentargets_summ$ensg))
length(exp_with_omim); print(length(unique(omim_genes)))

##### zscores
###SUMMARY STATS
eOutliers_blood_hg19<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Blood.hg19.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='blood',build='hg19')
eOutliers_fibroblast_hg19<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Fibroblast.hg19.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='fibroblast',build='hg19')
eOutliers_blood_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='blood',build='hg38')
eOutliers_fibroblast_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='fibroblast',build='hg38')
eOutliers_blood_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='blood',build='chm13')
eOutliers_fibroblast_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='fibroblast',build='chm13')
eOutliers_blood_chm13ensembl<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13ensembl.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='blood',build='chm13_ensembl')
eOutliers_fibroblast_chm13ensembl<-fread("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13ensembl.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.txt")%>%group_by(sample_id)%>%mutate(rank_ind=rank(desc(abs(zscore)),ties.method="first"),tissue='fibroblast',build='chm13_ensembl')

all_eoutliers=rbind(eOutliers_blood_hg19,eOutliers_fibroblast_hg19,
  eOutliers_blood_chm13,eOutliers_fibroblast_chm13,
  #eOutliers_blood_chm13ensembl,eOutliers_fibroblast_chm13ensembl,
  eOutliers_blood_hg38,eOutliers_fibroblast_hg38)
all_eoutliers$status<-sample_to_status[all_eoutliers$sample_id]
all_eoutliers$ensg<-gsub("\\..*","",all_eoutliers$gene)

myzscores=c("eOutliers_blood_hg19","eOutliers_fibroblast_hg19",
            "eOutliers_blood_hg38","eOutliers_fibroblast_hg38",
            "eOutliers_blood_chm13","eOutliers_fibroblast_chm13") #,
            #"eOutliers_blood_chm13ensembl","eOutliers_fibroblast_chm13ensembl")
summ_stats_exp_df=data.frame()
all_eoutliers=data.frame()

for(this_zscore in myzscores){
  my_build=sapply(str_split(this_zscore,"_"),"[[",3)
  for(over_under in c("Under eOutlier","Over eOutlier")){ #
    
    if(over_under=="Over eOutlier"){
      print(over_under)
      
      outliers=get(this_zscore)%>%dplyr::filter(abs(zscore)>3)
      all_eoutliers=rbind(all_eoutliers,cbind(outliers, exp_spl=over_under))
    }else{
      print(over_under)
      
      outliers=get(this_zscore)%>%dplyr::filter((zscore)<(-3))
      all_eoutliers=rbind(all_eoutliers,cbind(outliers, exp_spl=over_under))
    }
    
    num_outliers=outliers%>%nrow()
    num_outlier_genes=length(unique(outliers$gene ))
    uniq_omim_gene_overlap=outliers%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene %in% omim_genes)%>%pull(gene)%>%unique()%>%length()
    summ_stats_exp_df=rbind(as.data.frame(summ_stats_exp_df),
                            as.data.frame(cbind(this_zscore,num_outliers,num_outlier_genes,uniq_omim_gene_overlap, exp_spl=over_under,build=my_build)))
    }
  }

all_eoutliers$status<-sample_to_status[all_eoutliers$sample_id]

summstats_exp=summ_stats_exp_df%>%separate(this_zscore, c("tissue","build"),sep = "_")

###COMPARISON
## read files
eOutliers_blood_hg19_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz")%>% #%>%mutate(gene=gsub("\\..*","",gene))%>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%mutate(rank_build2=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(comparison="hg19:hg38",tissue="blood")
eOutliers_fibroblast_hg19_hg38<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg19.hg38.txt.gz")%>% #%>%mutate(gene=gsub("\\..*","",gene))%>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%mutate(rank_build2=rank(desc(abs(z_hg38)),ties.method="first"))%>%
  mutate(comparison="hg19:hg38",tissue="fibroblast")
eOutliers_blood_hg38_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg38.chm13ensembl.txt.gz")%>% #%>%mutate(gene=gsub("\\..*","",gene))%>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg38)),ties.method="first"))%>%mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg38:chm13",tissue="blood")
eOutliers_fibroblast_hg38_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg38.chm13ensembl.txt.gz")%>% #%>%mutate(gene=gsub("\\..*","",gene))%>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg38)),ties.method="first"))%>%mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg38:chm13",tissue="fibroblast")
eOutliers_blood_hg19_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Blood.dedupOptical_minMQ255.zscores_genesdiff.hg19.chm13ensembl.txt.gz")%>% #%>%mutate(gene=gsub("\\..*","",gene))%>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg19:chm13",tissue="blood")
eOutliers_fibroblast_hg19_chm13<-fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/Fibroblast.dedupOptical_minMQ255.zscores_genesdiff.hg19.chm13ensembl.txt.gz")%>% #%>%mutate(gene=gsub("\\..*","",gene))%>%
  group_by(sample_id)%>%mutate(rank_build1=rank(desc(abs(z_hg19)),ties.method="first"))%>%mutate(rank_build2=rank(desc(abs(z_chm13)),ties.method="first"))%>%
  mutate(comparison="hg19:chm13",tissue="fibroblast")




#combine so if different annotation we just note it but still compare
combine_diff_annotation<-function(z){
  nas_chm13=z%>%dplyr::filter(str_detect(gene,"-"))%>%mutate(gene2=sapply(strsplit(gene,"\\."),"[[",1))%>%dplyr::filter(is.na(z_diff))%>%pull(gene2)%>%unique()
  no_nas=z%>%dplyr::filter(!is.na(z_diff))%>%mutate(annotation_different=FALSE)
  nas<-z%>%dplyr::filter(is.na(z_diff))%>%rename(gene_withversion=gene)%>%
    mutate(gene=sapply(strsplit(gene_withversion,"\\."),"[[",1))%>%
    group_by(sample_id,gene)%>%
    mutate(n=n())
  keep_na=nas%>%dplyr::filter(n==1)%>%dplyr::select(-n, -gene_withversion)%>%mutate(annotation_different=NA)
  get_to_combine=nas%>%dplyr::filter(n==2)%>%dplyr::select(-n, -gene_withversion) %>%dplyr::filter(!(gene%in%nas_chm13))
  if(nrow(get_to_combine)>0){
    #combine_diff_versions=melt(get_to_combine)%>%dplyr::filter(!is.na(value))%>%group_by(sample_id,gene)%>%pivot_wider(names_from=variable,values_from=value)

    if("z_hg38"%in%colnames(get_to_combine))  {
      combine_diff_versions=get_to_combine%>%group_by(sample_id,gene,comparison,tissue)%>%
        summarise(z_hg19=max(z_hg19,na.rm = T),z_hg38=max(z_hg38,na.rm=T),rank_build1=max(rank_build1,na.rm=T),rank_build2=max(rank_build2,na.rm=T))
    }else{
      combine_diff_versions=get_to_combine%>%group_by(sample_id,gene,comparison,tissue)%>%
        summarise(z_hg38=max(z_hg38,na.rm = T),z_chm13=max(z_chm13,na.rm=T),rank_build1=max(rank_build1,na.rm=T),rank_build2=max(rank_build2,na.rm=T))
    }

    
    combine_diff_versions_z_diff=combine_diff_versions[,6]-combine_diff_versions[,5]
    combine_diff_versions$z_diff=combine_diff_versions_z_diff
    colnames(combine_diff_versions)[ncol(combine_diff_versions)]<-"z_diff"
    # combine_diff_versions_all=cbind(combine_diff_versions,combine_diff_versions_z_diff)
    combine_diff_versions=combine_diff_versions%>%mutate(annotation_different=TRUE)
    revamped_df=rbind(as.data.frame(no_nas),
                      as.data.frame(keep_na),
                      as.data.frame(combine_diff_versions))
  }else{
    print("not any needed to recombine")
    revamped_df=rbind(as.data.frame(no_nas),
                      as.data.frame(keep_na))
  }
  


  return(revamped_df)
}
eOutliers_fibroblast_hg19_hg38_anno=combine_diff_annotation(eOutliers_fibroblast_hg19_hg38)
eOutliers_fibroblast_hg38_chm13_anno=combine_diff_annotation(eOutliers_fibroblast_hg38_chm13)
# eOutliers_fibroblast_hg19_chm13_anno=combine_diff_annotation(eOutliers_fibroblast_hg19_chm13)
eOutliers_blood_hg19_hg38_anno=combine_diff_annotation(eOutliers_blood_hg19_hg38)
eOutliers_blood_hg38_chm13_anno=combine_diff_annotation(eOutliers_blood_hg38_chm13)
# eOutliers_blood_hg19_chm13_anno=combine_diff_annotation(eOutliers_blood_hg19_chm13)

#outliers in both/outliers in hg19
get_outlier_prop<-function(z, build_reference,build_to_compare,tissue,omim_genes,chrY_chm13_genes,outlier_threshold=3){
  z=z%>%mutate(gene=gsub("\\..*","",gene))
  if(build_reference=="z_chm13" | build_to_compare=="z_chm13"){
    print("chm13")
    print(nrow(z))
    z=z%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))
    print(nrow(z))
  }
  print(paste0("What proportion of ",build_reference, " outliers are consistent in ",build_to_compare,"?"))
  num_outliers_consistent=round(z%>%dplyr::filter(get(build_reference)>=3 & get(build_to_compare)>=3)%>%nrow()/
                                  z%>%dplyr::filter(get(build_reference)>=3)%>%nrow()*100,3)
  print(paste0("percent of outliers that are consistent: ",num_outliers_consistent,"%"))
  num_outliers_notseenother=round(z%>%dplyr::filter(get(build_reference)>=3 & is.na(get(build_to_compare)))%>%nrow()/
                                    z%>%dplyr::filter(get(build_reference)>=3)%>%nrow()*100,3)
  print(paste0("percent not seen in the other build: ",num_outliers_notseenother ,"%"))
  for_df<-c(paste0(build_reference,"-",build_to_compare),tissue,
            num_outliers_consistent,
            num_outliers_notseenother)
  for(i in c(1,2,3)){
    denom=z%>%dplyr::filter(get(build_reference)>=outlier_threshold)%>%nrow()
    filt_zs=z%>%dplyr::filter(get(build_reference)>=outlier_threshold & get(build_to_compare)<i)
    # in_omim=filt_zs%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene%in%omim_genes)%>%pull(gene)%>%unique()%>%length()
    numerator=filt_zs%>%nrow()
    frac=numerator/denom*100
    print(paste0("z<",i, " is ", round(frac,outlier_threshold), "%"))
    # for_df<-c(for_df,round(frac,3),in_omim)
    for_df<-c(for_df,round(frac,3))
  }
  
  return_df=as.data.frame(t(for_df))
  colnames(return_df)<-c("comparison","tissue","consistent","not_in_other_build","z<1",
                         "z<2","z<3")
  # colnames(return_df)<-c("comparison","tissue","consistent","not_in_other_build","z<1", "omim_z<1",
  #                        "z<2", "omim_z<2","z<3", "omim_z<3")
  return(return_df)
}
get_outlier_consistency_df<-function(){
  eOutlier_consistency_df=data.frame()
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_blood_hg19_hg38_anno,"z_hg19","z_hg38","blood",omim_genes,chrY_chm13_genes))
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_blood_hg19_hg38_anno,"z_hg38","z_hg19","blood",omim_genes,chrY_chm13_genes))
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_fibroblast_hg19_hg38_anno,"z_hg19","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_fibroblast_hg19_hg38_anno,"z_hg38","z_hg19","fibroblast",omim_genes,chrY_chm13_genes))
  
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_blood_hg38_chm13_anno,"z_hg38","z_chm13","blood",omim_genes,chrY_chm13_genes))
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_blood_hg38_chm13_anno,"z_chm13","z_hg38","blood",omim_genes,chrY_chm13_genes))
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_fibroblast_hg38_chm13_anno,"z_hg38","z_chm13","fibroblast",omim_genes,chrY_chm13_genes))
  eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(eOutliers_fibroblast_hg38_chm13_anno,"z_chm13","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
  
  # eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(blood_hg19_chm13_anno,"z_hg19","z_chm13","blood",omim_genes,chrY_chm13_genes))
  # eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(blood_hg19_chm13_anno,"z_chm13","z_hg19","blood",omim_genes,chrY_chm13_genes))
  # eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(fibroblast_hg19_chm13_anno,"z_hg19","z_chm13","fibroblast",omim_genes,chrY_chm13_genes))
  # eOutlier_consistency_df=rbind(eOutlier_consistency_df,get_outlier_prop(fibroblast_hg19_chm13_anno,"z_chm13","z_hg19","fibroblast",omim_genes,chrY_chm13_genes))
  return(eOutlier_consistency_df)
} 
eOutlier_consistency_df=get_outlier_consistency_df()

make_zdiff_plot<-function(eOutlier_consistency_df,sOutlier_consistency_df){
  outlier_consistency_df=rbind(cbind(sOutlier_consistency_df,assay='splicing'),
                               cbind(eOutlier_consistency_df,assay='expression'))
  
  Outlier_consistency_df_reformat=outlier_consistency_df%>%mutate(comparison=gsub('z_','',comparison))%>%mutate(comparison=gsub('-',':',comparison))%>%
    separate(comparison,c('reference_build','comparison_build'))%>%
    dplyr::filter((reference_build=='hg38'&comparison_build=='hg19')|(reference_build=='chm13'&comparison_build=='hg38'))%>%
    mutate(comparison=paste0(comparison_build,':',reference_build))%>%
    #mutate('original_z<3'=as.numeric(as.character(get('z<3'))))%>%
    mutate('z<2'=as.numeric(as.character(get('z<2')))-as.numeric(as.character(get('z<1'))),
           'z<3'=as.numeric(as.character(get('z<3')))-as.numeric(as.character(get('z<2')))-as.numeric(as.character(get('z<1'))))
  Outlier_consistency_df_reformat_melted=melt(Outlier_consistency_df_reformat,id=c('reference_build','comparison_build','tissue','comparison','assay'))%>%
    dplyr::filter(variable!='omim_z<1'&variable!='omim_z<2'&variable!='omim_z<3')%>%mutate(value=as.numeric(value)) %>%
    mutate(variable=ifelse(variable!='not_in_other_build',as.character(variable),'not in other build'))
  Outlier_consistency_df_reformat_melted$variable=factor(Outlier_consistency_df_reformat_melted$variable,levels=c('consistent','z<3','z<2','z<1','not in other build'))
  Outlier_consistency_df_reformat_melted$variable=factor(Outlier_consistency_df_reformat_melted$variable,levels=c('consistent','z<3','z<2','z<1','not in other build'))
  ggplot(Outlier_consistency_df_reformat_melted,aes(x=value,y=comparison,fill=comparison,alpha=variable))+
    scale_fill_manual(values=c("#d07c10","#006B2F"))+
    # scale_alpha_discrete(range=c(0.4,1))+
    scale_alpha_manual(values = c(1,0.85,0.65,0.45,0.3,0))+
    
    # coord_flip() + scale_x_continuous(trans='log10')+
    #scale_x_log10(breaks=c(0.01,0.1,1),labels=c(0.01,0.1,1))+
    #scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+   # xlim(values=c(0,101))+
    geom_bar(stat = 'identity',color="#363738")+
    #xlim(values=c(0,1.1))+
    # scale_fill_manual(values=c("#37771D","#609c49","#92C47C","#D9EAD3","#9E9E9E"))+
    theme_classic(base_size=15)+
    facet_wrap(~assay*tissue)
  ggplot(Outlier_consistency_df_reformat_melted%>%dplyr::filter(variable!='consistent'),aes(x=value,y=comparison,fill=comparison,alpha=variable))+
    scale_fill_manual(values=c("#d07c10","#006B2F"))+
    # scale_alpha_discrete(range=c(0.4,1))+
    scale_alpha_manual(values = c(0.85,0.65,0.45,0.3,0))+
    
    # coord_flip() + scale_x_continuous(trans='log10')+
    #scale_x_log10(breaks=c(0.01,0.1,1),labels=c(0.01,0.1,1))+
    #scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+   # xlim(values=c(0,101))+
    geom_bar(stat = 'identity',color="#363738")+
    #xlim(values=c(0,1.1))+
    # scale_fill_manual(values=c("#37771D","#609c49","#92C47C","#D9EAD3","#white"))+
    theme_classic(base_size=15)+
    facet_wrap(~assay*tissue)
  
  
}
make_zdiff_plot(eOutlier_consistency_df,sOutlier_consistency_df)

  #######pull the OMIM genes w large zdiff (z>3, z<1)
get_omim_list_exp<-function(z,outlier_threshold=3,build_reference,build_to_compare,tissue,omim_genes,chrY_chm13_genes){
    if(build_reference=="z_chm13" | build_to_compare=="z_chm13"){
      print("chm13")
      print(nrow(z))
      z=z%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(!(gene %in%chrY_chm13_genes$gene))
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
get_eOutlier_omim_diff_df<-function(){
  eOutlier_omim_diff_df=data.frame()
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(blood_hg38_chm13,3,"z_hg38","z_chm13","blood",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(fibroblast_hg38_chm13,3,"z_hg38","z_chm13","fibroblast",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(blood_hg19_hg38,3,"z_hg19","z_hg38","blood",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(fibroblast_hg19_hg38,3,"z_hg19","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(blood_hg38_chm13,3,"z_chm13","z_hg38","blood",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(fibroblast_hg38_chm13,3,"z_chm13","z_hg38","fibroblast",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(blood_hg19_hg38,3,"z_hg38","z_hg19","blood",omim_genes,chrY_chm13_genes))
  eOutlier_omim_diff_df=rbind(eOutlier_omim_diff_df,get_omim_list_exp(fibroblast_hg19_hg38,3,"z_hg38","z_hg19","fibroblast",omim_genes,chrY_chm13_genes))
  return(eOutlier_omim_diff_df)
}
eOutlier_omim_diff_df=get_eOutlier_omim_diff_df()



#overlap w solved
solved_and_candidates<-fread("/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_AND_candidate_dec2022.txt")
exp_omim_diff_solvedgene=eOutlier_omim_diff_df%>%dplyr::filter(gene%in%solved_and_candidates$ENSG)
blood_hg38_chm13%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene%in%exp_omim_diff_solvedgene$gene)%>%dplyr::filter(abs(z_diff)>2)%>%mutate(status=sample_to_status[sample_id])



#top genes ranks comparison
all_comparisons=rbind(blood_hg19_hg38,fibroblast_hg19_hg38,
                      blood_hg38_chm13,fibroblast_hg38_chm13) #,
                     # blood_hg19_chm13,fibroblast_hg19_chm13)
all_comparisons_top20=all_comparisons%>%dplyr::filter(rank_build1<20 | rank_build2<20)

ggplot(all_comparisons_top20,
       aes(x=rank_build1,y=rank_build2,alpha=0.3,color=tissue))+
  geom_point()+
  scale_x_continuous(trans="log10")+
  theme_bw()+
  stat_regline_equation(label.x=3,aes(label = ..rr.label..)) +
  ggtitle("Expression Outlier Rank Comparison")+
  facet_wrap(~comparison)

## diff exp overlap
mytissues=c("Blood","Fibroblast","Muscle","PBMC","iPSC_NPC","iPS")
all_diffexp=data.frame()
all_diffexp_anylogfc=data.frame()
for(this_tissue in mytissues){
  for(this_buildcomp in c("chm13ensembl_vs_hg38","chm13ensembl_vs_hg19","hg38_vs_hg19")){
    my_filename=paste0("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/",this_tissue,".dedupOptical_minMQ255.diff_expression_",this_buildcomp,".limma_voom_dream.significant.txt.gz")
    mydf<-fread(my_filename)%>%dplyr::filter(abs(logFC)>1)%>%mutate(tissue=this_tissue, comparison=this_buildcomp)
    all_diffexp=rbind(all_diffexp,mydf)
    
    my_filename_anylogfc=paste0("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/",this_tissue,".dedupOptical_minMQ255.diff_expression_",this_buildcomp,".limma_voom_dream.all.txt.gz")
    mydf_anylogfc<-fread(my_filename_anylogfc)%>%mutate(tissue=this_tissue, comparison=this_buildcomp)
    all_diffexp_anylogfc=rbind(all_diffexp_anylogfc,mydf_anylogfc)
    
  }
  
}


##opentargets
all_diffexp_anno<-merge(all_diffexp,opentargets_summ,by.x="gene",by.y="ensg",all.x = T)
all_diffexp_anno<-merge(all_diffexp_anno,unique(omim_df[,c("ensg",'confidence')]),by.x="gene",by.y="ensg",all.x = T)
ggplot(all_diffexp_anno,aes(x=logFC,y=target_score,color=tissue,alpha=0.8))+
  geom_point()+
  scale_color_manual(values=c("#D90025","#27A795","#0470B5","#BD71DC","#AACA2F","#EF7C18"))+
  facet_wrap(~comparison)+
  theme_bw()
COL1A1=fibroblast_hg38_chm13%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene=="ENSG00000108821")
KRAS=fibroblast_hg38_chm13%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene=="ENSG00000133703")
RB1=fibroblast_hg38_chm13%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene=="ENSG00000139687")

compare_outliers_diffexp<-function(z,comp1,comp2,diff_exp,omim_genes,my_tissue,my_comparison){
  diff_exp=diff_exp%>%dplyr::filter(tissue==my_tissue & comparison==my_comparison)%>%dplyr::filter(abs(logFC)>1 & adj.P.Val<0.05)
  z_diff_exp<-z%>%mutate(diff_exp=ifelse(gene %in% diff_exp$gene,"differentially_quantified","not_differentially_quantified"))
  z_diff_exp_outliers<-z_diff_exp%>%dplyr::filter(abs(get(comp1))>3 | abs(get(comp2))>3)
  z_diff_exp_largediff<-z_diff_exp%>%dplyr::filter((abs(get(comp1))>3 & abs(get(comp2))<1 ) |
                                                    (abs(get(comp1))<1 & abs(get(comp2))>3 )
                                                   )
  uniq_genes_quantorno=z_diff_exp_largediff%>%ungroup()%>%dplyr::select(gene,diff_exp)%>%unique()
  mytable=table(uniq_genes_quantorno$diff_exp)
  mydf=cbind(comparison=my_comparison,tissue=my_tissue,
             total_diffexp_genes=nrow(diff_exp),
        diffexp_genes_no_outliers=table(diff_exp$gene %in% z_diff_exp_outliers$gene)["FALSE"],
        diffexp_gene_yes_outliers=table(diff_exp$gene %in% z_diff_exp_outliers$gene)["TRUE"],
        uniq_gene_outliers=length(unique(z_diff_exp_outliers$gene)),
        num_large_differences=length(unique(z_diff_exp_largediff$gene)),
        t(data.frame(unclass(mytable))))
  
  # uniq_omim_gene_overlap=outliers%>%mutate(gene=gsub("\\..*","",gene))%>%dplyr::filter(gene %in% omim_genes)%>%pull(gene)%>%unique()%>%length()
  return(as.data.frame(mydf))
}

outlier_diffexp_table=rbind(compare_outliers_diffexp(eOutliers_blood_hg19_hg38,"z_hg19","z_hg38",all_diffexp,omim_genes,"Blood","hg38_vs_hg19"),
                            compare_outliers_diffexp(eOutliers_fibroblast_hg19_hg38,"z_hg19","z_hg38",all_diffexp,omim_genes,"Fibroblast","hg38_vs_hg19"),
                            compare_outliers_diffexp(eOutliers_blood_hg38_chm13,"z_hg38","z_chm13",all_diffexp,omim_genes,"Blood","chm13ensembl_vs_hg38"),
                            compare_outliers_diffexp(eOutliers_fibroblast_hg38_chm13,"z_hg38","z_chm13",all_diffexp,omim_genes,"Fibroblast","chm13ensembl_vs_hg38"))
all_outliers_comparisons=rbind(eOutliers_blood_hg19_hg38,eOutliers_fibroblast_hg19_hg38,eOutliers_blood_hg38_chm13,eOutliers_fibroblast_hg38_chm13)
combine_outliers_diff_exp<-function(zcomp,diffexp){
  comparison_dic=c("hg38_vs_hg19","chm13ensembl_vs_hg38")
  names(comparison_dic)=c("hg19:hg38","hg38:chm13")
  
  zcomp_summ=zcomp%>%group_by(gene,comparison,tissue)%>%summarise(max_zdiff=max(abs(z_diff)),ave_abs_zdiff=mean(abs(z_diff)))%>%
    mutate(gene=gsub("\\..*","",gene))%>%
    mutate(tissue=str_to_title(tissue))%>%
    mutate(comparison_notation=comparison,comparison=comparison_dic[comparison]) #%>%head()
  z_and_diffexp=merge(zcomp_summ,diffexp,by=c("gene","comparison","tissue"))
  return(z_and_diffexp)
}
zsumm_and_diffexp=combine_outliers_diff_exp(all_outliers_comparisons,all_diffexp_anylogfc)
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

make_zsumm_diffexp_plot<-function(zsumm_and_diffexp){
  library(ggpubr)
  zsumm_and_diffexp=zsumm_and_diffexp%>%mutate(abs_logFC=abs(logFC))
  ggscatter(data=(zsumm_and_diffexp),x="abs_logFC",y="max_zdiff",color="adj.P.Val",alpha=0.4, add = "reg.line")+
    theme_bw()+
    facet_wrap(~comparison_notation*tissue,scales = "free_y")+
    stat_cor(
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
      label.x = 5.5,label.y=0.25
    )+
  # stat_cor(label.y = .25,label.x=5.5)+
    scale_color_gradient2(high="#ba54b7",low="#37827b",mid="#06b8c2",midpoint=0.01)
  return(my_plot)
}



get_solved_outlier<-function(){
  # compare_outliers_diffexp(blood_hg38_chm13,"z_hg38","z_chm13",diffexp_blood_hg38_chm13,omim_genes,"blood","hg38:chm13"),## Solved cases
  solved_stanford_transcriptomics<-fread("/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases_blood_fibroblast.txt") 
  rdid_to_ensg<-solved_stanford_transcriptomics$ENSG
  names(rdid_to_ensg)<-solved_stanford_transcriptomics$RDID
  solved_stanford_transcriptomics_blood<-solved_stanford_transcriptomics%>%dplyr::filter(Tissue=="Blood")
  solved_stanford_transcriptomics_fibroblast<-solved_stanford_transcriptomics%>%dplyr::filter(Tissue=="Fibroblast")
  #dplyr::filter(z_diff_fibroblast_sorted,z_hg38!=0 & z_hg19 !=0)
  get_solved_genes<-function(df,rdid_to_ensg){
    return(df%>%mutate(ensg=str_replace(gene,"\\..*$",""))%>%dplyr::filter(sample_id %in% names(rdid_to_ensg) & rdid_to_ensg[sample_id]==ensg))
  }
  
  # stanford_solved_fibroblast_hg19_chm13<-get_solved_genes(fibroblast_hg19_chm13_anno,rdid_to_ensg)
  stanford_solved_fibroblast_hg38_chm13<-get_solved_genes(eOutliers_fibroblast_hg38_chm13,rdid_to_ensg)
  stanford_solved_fibroblast_hg19_hg38<-get_solved_genes(eOutliers_fibroblast_hg19_hg38,rdid_to_ensg)
  # stanford_solved_blood_hg19_chm13<-get_solved_genes(blood_hg19_chm13_anno,rdid_to_ensg)
  stanford_solved_blood_hg38_chm13<-get_solved_genes(eOutliers_blood_hg38_chm13,rdid_to_ensg)
  stanford_solved_blood_hg19_hg38<-get_solved_genes(eOutliers_blood_hg19_hg38,rdid_to_ensg)
  
  stanford_solved_hg19_hg38=rbind(stanford_solved_blood_hg19_hg38%>%mutate(tissue="Blood"),
                                  stanford_solved_fibroblast_hg19_hg38%>%mutate(tissue="Fibroblast"))
  # stanford_solved_hg19_chm13=rbind(stanford_solved_blood_hg19_chm13%>%mutate(tissue="Blood"),
  #                                  stanford_solved_fibroblast_hg19_chm13%>%mutate(tissue="Fibroblast"))
  stanford_solved_hg38_chm13=rbind(stanford_solved_blood_hg38_chm13%>%mutate(tissue="Blood"),
                                   stanford_solved_fibroblast_hg38_chm13%>%mutate(tissue="Fibroblast"))
  return(list(stanford_solved_hg19_hg38,stanford_solved_hg38_chm13))
}
get_solved_outliers=get_solved_outlier()
toplot=rbind(get_solved_outliers[[2]],get_solved_outliers[[1]])%>%mutate(z_comp2=ifelse(is.na(z_chm13),z_hg19,z_chm13))
myplot=ggplot(toplot,aes(x=z_comp2,y=z_hg38,alpha=0.8,color=abs(z_diff),shape=tissue))+
  geom_abline(intercept=0,slope=1,color="red")+
  scale_color_gradient(low="#5de3ab",high="#280a47")+
  facet_wrap(~comparison)+
  ggtitle("Solved Stanford Z-Scores ")+theme_bw(base_size = 18)+xlim(c(-4,4))+ylim(c(-4,4))+
  geom_point(aes(alpha=0.8),size=3) #+facet_grid(~tissue,scales="free_x")
myplot

solved_outlier_diffexp<-function(){
  solved_diffexp_blood_hg38_chm13=diffexp_blood_hg38_chm13[diffexp_blood_hg38_chm13$gene%in%solved_stanford_transcriptomics$ENSG,]
  solved_diffexp_blood_hg19_hg38=diffexp_blood_hg19_hg38[diffexp_blood_hg19_hg38$gene%in%solved_stanford_transcriptomics$ENSG,]
  solved_diffexp_blood_hg19_chm13=diffexp_blood_hg19_chm13[diffexp_blood_hg19_chm13$gene%in%solved_stanford_transcriptomics$ENSG,]
  solved_diffexp_fibroblast_hg38_chm13=diffexp_fibroblast_hg38_chm13[diffexp_fibroblast_hg38_chm13$gene%in%solved_stanford_transcriptomics$ENSG,]
  solved_diffexp_fibroblast_hg19_hg38=diffexp_fibroblast_hg19_hg38[diffexp_fibroblast_hg19_hg38$gene%in%solved_stanford_transcriptomics$ENSG,]
  solved_diffexp_fibroblast_hg19_chm13=diffexp_fibroblast_hg19_chm13[diffexp_fibroblast_hg19_chm13$gene%in%solved_stanford_transcriptomics$ENSG,]
  #omim diff exp
  omim_diffexp_blood_hg38_chm13=diffexp_blood_hg38_chm13[sapply(strsplit(diffexp_blood_hg38_chm13$gene,"\\."),"[[",1)%in%omim_genes,]; print(nrow(omim_diffexp_blood_hg38_chm13))
  omim_diffexp_blood_hg19_chm13=diffexp_blood_hg19_chm13[sapply(strsplit(diffexp_blood_hg19_chm13$gene,"\\."),"[[",1)%in%omim_genes,]; print(nrow(omim_diffexp_blood_hg19_chm13))
  omim_diffexp_blood_hg19_hg38=diffexp_blood_hg19_hg38[sapply(strsplit(diffexp_blood_hg19_hg38$gene,"\\."),"[[",1)%in%omim_genes,]; print(nrow(omim_diffexp_blood_hg19_hg38))
  omim_diffexp_fibroblast_hg38_chm13=diffexp_fibroblast_hg38_chm13[sapply(strsplit(diffexp_fibroblast_hg38_chm13$gene,"\\."),"[[",1)%in%omim_genes,]; print(nrow(omim_diffexp_fibroblast_hg38_chm13))
  omim_diffexp_fibroblast_hg19_chm13=diffexp_fibroblast_hg19_chm13[sapply(strsplit(diffexp_fibroblast_hg19_chm13$gene,"\\."),"[[",1)%in%omim_genes,]; print(nrow(omim_diffexp_fibroblast_hg19_chm13))
  omim_diffexp_fibroblast_hg19_hg38=diffexp_fibroblast_hg19_hg38[sapply(strsplit(diffexp_fibroblast_hg19_hg38$gene,"\\."),"[[",1)%in%omim_genes,]; print(nrow(omim_diffexp_fibroblast_hg19_hg38))
  
  
  # blood_hg19_hg38_omim_diffexp=blood_hg19_hg38%>%mutate(gene2=sapply(strsplit(gene,"\\."),"[[",1))%>%dplyr::filter(gene2 %in% omim_diffexp_blood_hg19_hg38$gene )
  blood_hg19_hg38_omim_diffexp=blood_hg19_hg38%>%dplyr::filter(gene %in% omim_diffexp_blood_hg19_hg38$gene )
  fibroblast_hg19_hg38_omim_diffexp=fibroblast_hg19_hg38%>%dplyr::filter(gene %in% omim_diffexp_fibroblast_hg19_hg38$gene )
  blood_hg38_chm13_omim_diffexp=blood_hg38_chm13%>%dplyr::filter(gene %in% omim_diffexp_blood_hg38_chm13$gene )
  fibroblast_hg38_chm13_omim_diffexp=fibroblast_hg38_chm13%>%dplyr::filter(gene %in% omim_diffexp_fibroblast_hg38_chm13$gene )
  
  fibroblast_hg19_hg38_omim_diffexp%>%dplyr::filter(!is.na(z_diff))%>%arrange(desc(abs(z_diff)))%>%pull(gene2)
  a=blood_hg38_chm13_omim_diffexp%>%dplyr::filter(!is.na(z_diff))%>%arrange(desc(abs(z_diff))) #%>%pull(gene) #%>%head(15)%>%tail(10)
  chm13_s13genes=fread("/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/s13_medicallyrelevant_chm13_impacted.csv",header=T) #%>%pull(V1)
  chm13genesdiffs=chm13_s13genes$gene_name[chm13_s13genes$gene_name%in%ensg_to_hgnc[sapply(strsplit(a$gene,"\\."),"[[",1)]]
  chm13_s13genes%>%dplyr::filter(gene_name%in%chm13genesdiffs&`Impacted?`>0)%>%pull(gene_name)
  b=diffexp_blood_hg38_chm13%>%dplyr::mutate(gene2=sapply(strsplit(gene,"\\."),"[[",1))%>%dplyr::filter(gene2 %in% hgnc_to_ensg[chm13_s13genes%>%dplyr::filter(gene_name%in%chm13genesdiffs&`Impacted?`>0)%>%pull(gene_name)]) %>%arrange(desc(abs(logFC)))
  
  
  cdh23=blood_hg38_chm13_omim_diffexp%>%filter(gene=="ENSG00000107736.21")
  cdh23_melted<-cdh23%>%pivot_longer(cols = z_hg38:z_chm13)
  cdh23_melted$name<-factor(cdh23_melted$name,levels=c("z_hg38","z_chm13"))
  samples_that_changed=cdh23%>%dplyr::filter((abs(z_hg38)>3&abs(z_chm13)<3) | (abs(z_hg38)<3&abs(z_chm13)>3)) %>%pull(sample_id)
  ggplot(cdh23_melted,aes(x=name,y=value))+theme_bw()+
    scale_colour_gradient2(mid="#e8df76",low="#e8df76",high="#E98301")+ylab("z-score")+
    geom_line(aes(group=sample_id,color=abs(z_diff),alpha=0.5))+
    #geom_line(data=cdh23_melted%>%dplyr::filter(sample_id%in%samples_that_changed),
    #          aes(group=sample_id))+
    geom_point(aes(alpha=0.5))
  
}
solved_outlier_diffexp_plot=get_solved_outlier_diffexp()


##oberlap w main table
main_table=fread('/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/all_genes_annotated_stupid_big.diffexp.mappingstats.outliers.updated_Oct2023.tsv')
main_table=main_table%>%mutate(hg38_chm13_annotation_comparison_edit=ifelse(str_detect(hg38_chm13_annotation_comparison,"hg38"),"hg38_spec_gene",
                                                                            ifelse(str_detect(hg38_chm13_annotation_comparison,"chm13"),'chm13_spec_gene',hg38_chm13_annotation_comparison)))%>%
  mutate(any_clin_relevant=ifelse(is_opentargets_gene==1|is_omim_gene==1|is_cancer_gene==1|is_clinvar_gene==1, 
                                  "clinically_relevant","no clinical relevance"))
hg19_spec_genes=main_table%>%dplyr::filter(hg19_hg38_annotation_comparison=='hg19_specific')%>%pull(gene_id)%>%unique()
hg38_compared2hg19_spec_genes=main_table%>%dplyr::filter(hg19_hg38_annotation_comparison=='hg38_specific')%>%pull(gene_id)%>%unique()
hg38_compared2chm13_spec_genes=main_table%>%dplyr::filter(hg38_chm13_annotation_comparison_edit=='hg38_spec_gene')%>%pull(gene_id)%>%unique()
chm13_spec_genes=main_table%>%dplyr::filter(hg38_chm13_annotation_comparison_edit=='chm13_spec_gene')%>%pull(gene_id)%>%unique()
zthresh=3
hg19_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%hg19_spec_genes&abs(zscore)>zthresh)%>%mutate(unique_in="hg19")
hg38compared2hg19_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%hg38_compared2hg19_spec_genes&abs(zscore)>zthresh)%>%mutate(unique_in="hg38_ascomparedtohg19")%>%dplyr::filter(build=='hg38')
hg38compared2chm13_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%hg38_compared2chm13_spec_genes&abs(zscore)>zthresh)%>%mutate(unique_in="hg38_ascomparedtochm13")%>%dplyr::filter(build=='hg38')
chm13_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%chm13_spec_genes&abs(zscore)>zthresh)%>%mutate(unique_in="chm13")

hg19_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%hg19_spec_genes&rank_ind<=20&status=="Case")%>%mutate(unique_in="hg19")%>%mutate(comparison="hg19:hg38")
hg38compared2hg19_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%hg38_compared2hg19_spec_genes&rank_ind<=20&status=="Case")%>%mutate(unique_in="hg38_ascomparedtohg19")%>%dplyr::filter(build=='hg38')%>%mutate(comparison="hg19:hg38")
hg38compared2chm13_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%hg38_compared2chm13_spec_genes&rank_ind<=20&status=="Case")%>%mutate(unique_in="hg38_ascomparedtochm13")%>%dplyr::filter(build=='hg38')%>%mutate(comparison="hg38:chm13")
chm13_spec_eoutliers=all_eoutliers%>%dplyr::filter(ensg%in%chm13_spec_genes&rank_ind<=20&status=="Case")%>%mutate(unique_in="chm13")%>%mutate(comparison="hg38:chm13")

all_spec_eoutliers_top20=rbind(hg19_spec_eoutliers,hg38compared2hg19_spec_eoutliers,hg38compared2chm13_spec_eoutliers,chm13_spec_eoutliers)
all_nonissues=get_anno_spec_eoutlier_info(all_spec_eoutliers_top20,'hg19',main_table)

get_anno_spec_eoutlier_info<-function(myoutliers,build,main_table){
   # myoutliers=myoutliers%>%dplyr::filter(rank_ind<=20)
  uniq_outliers=unique(myoutliers%>%pull(gene))
  num_outliers=length(uniq_outliers)
  num_outliers_control=length(unique(myoutliers%>%dplyr::filter(status=="Control")%>%pull(gene)))
  outliers_case=myoutliers%>%dplyr::filter(status=="Case")
  outliers_control=myoutliers%>%dplyr::filter(status=="Control")
  
  print('clinically relevant?')
  print((main_table%>%dplyr::filter(gene_id %in% myoutliers$ensg & any_clin_relevant=='clinically_relevant')%>%dplyr::select(gene_id,gene_name)))
  print(table(main_table%>%dplyr::filter(gene_id %in% myoutliers$ensg )%>%dplyr::select(tissue,any_clin_relevant)))

  num_outliers_case=length(unique(outliers_case$gene))
  print(paste0("num outliers: ", num_outliers," num outliers case only: ",num_outliers_case, " control only: ",num_outliers_control))
  print('case')
  genetable_outliers_cases=main_table%>%dplyr::filter(gene_id%in%outliers_case$ensg)
  print(genetable_outliers_cases%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
    dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table())
  nonissue_genes=genetable_outliers_cases%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
    dplyr::filter(any_issue=='non_issue')
  print('control')
  genetable_outliers_controls=main_table%>%dplyr::filter(gene_id%in%outliers_control$ensg)
  print(genetable_outliers_controls%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
          dplyr::select(any_issue,gene_id)%>%unique()%>%pull(any_issue)%>%table())
  nonissue_genes=genetable_outliers_controls%>%mutate(any_issue=ifelse(!is.na(get(paste0(build,"_issue")))|!is.na(get(paste0(build,"_blacklist_reason"))),"issue","non_issue"))%>%
    dplyr::filter(any_issue=='non_issue')
  return(nonissue_genes)
  
}
hg19_spec_nonissue=get_anno_spec_eoutlier_info(hg19_spec_eoutliers,'hg19',main_table)
hg38compared2hg19_spec_nonissue=get_anno_spec_eoutlier_info(hg38compared2hg19_spec_eoutliers,'hg38',main_table)
hg38compared2chm13_spec_nonissue=get_anno_spec_eoutlier_info(hg38compared2chm13_spec_eoutliers,'hg38',main_table)
chm13_spec_nonissue=get_anno_spec_eoutlier_info(chm13_spec_eoutliers,'chm13',main_table)
spec_eoutliers=rbind(hg19_spec_eoutliers,hg38compared2hg19_spec_eoutliers,hg38compared2chm13_spec_eoutliers,chm13_spec_eoutliers)
library(ggridges)

ggplot(spec_eoutliers,aes(x=rank_ind,fill=build,y=unique_in,alpha=0.3))+
  scale_fill_manual(values=c("#8C1515","#dea700","#01427a"))+
  scale_x_continuous(trans='log10')+
  
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
                      limits = c(1,7000),trim=TRUE)+
  theme_bw()+
  geom_vline(xintercept = 20)

## cases
solved_cases=fread("/oak/stanford/groups/smontgom/tannerj/UDN_ReferenceComparison_T2T/solved_cases/Solved_cases.phenotype_derived_genelist.transcriptome_prioritized.hg19_vs_hg38_vs_CHM13.txt")
solved_cases_ranks=solved_cases%>%dplyr::filter(rank_all<(50) |(SOLVED_GENE==TRUE&rank_all<100))
sig_dif_gene=main_table%>%dplyr::filter(adj_pval_chm13_hg38<0.05&tissue=="Blood"&abs(logFC_chm13_hg38)>0)%>%pull(gene_id)
to_plot=solved_cases_ranks%>%dplyr::filter(UDNID=="XX"&Tissue=="Blood")%>%dplyr::select(-zscore,-rank_underexp,-rank_overexp)%>%pivot_wider(names_from=BUILD,values_from=rank_all)%>%
  mutate(HG19=ifelse(is.na(HG19),-1,HG19),HG38=ifelse(is.na(HG38),-1,HG38),CHM13=ifelse(is.na(CHM13),-1,CHM13))%>%
  mutate(is_sig_dif_quant_hg38_chm13=ifelse(gene%in%sig_dif_gene,"significantly differentially quantified hg38:chm13","NOT significantly differentially quantified hg38:chm13 "))%>%
  mutate(quantification=ifelse(gene%in%hg38_compared2chm13_spec_genes,"hg38 annotation-specific",
                               ifelse(gene%in%chm13_spec_genes,"chm13 annotation-specific",
                                      ifelse(gene%in%sig_dif_gene, "significantly differentially quantified","no differential quantification")
                                      )
                               ))
  # mutate(quantification=ifelse(gene%in%hg38_compared2hg19_spec_genes,"hg19 annotation-specific",
  #                              ifelse(gene%in%hg19_spec_genes,"hg19 annotation-specific",
  #                                     ifelse(gene%in%sig_dif_gene, "significantly differentially quantified","no differential quantification")
  #                              )
  # ))
ggplot(to_plot,aes(x=HG38,y=CHM13,color=quantification,shape=SOLVED_GENE,label=GENE_SYMBOL))+
  theme_bw(base_size = 12)+
  facet_wrap(~Tissue,scales="free")+
  geom_text(hjust=0.5, vjust=-0.7)+
  scale_color_manual(values=c("#dea700","#01427a","black","#006B2F"))+
  scale_color_manual(values=c("black","#006B2F"))+
  geom_point(size=3)+ggtitle("Top 50 Genes: NEK9 (hg38:chm13)")


