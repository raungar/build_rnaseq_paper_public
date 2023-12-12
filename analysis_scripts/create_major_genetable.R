library(data.table)
library(tidyverse)

get_diffexp<-function(){
  ## diff exp overlap
  mytissues=c("Blood","Fibroblast","Muscle","PBMC","iPSC_NPC","iPS")
  all_diffexp=data.frame()
  all_diffexp_anylogfc=data.frame()
  for(this_tissue in mytissues){
    for(this_buildcomp in c("chm13ensembl_vs_hg38","chm13ensembl_vs_hg19","hg38_vs_hg19")){
      my_filename=paste0("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/",this_tissue,".dedupOptical_minMQ255.diff_expression_",this_buildcomp,".limma_voom_dream.significant.txt.gz")
      mydf<-fread(my_filename)%>%dplyr::filter(abs(logFC)>1)%>%mutate(tissue=this_tissue, comparison=this_buildcomp)
      all_diffexp=rbind(all_diffexp,mydf)
      
      my_filename_anylogfc=paste0("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/03_LIMMA_DiffExp_20221113/",this_tissue,".dedupOptical_minMQ255.diff_expression_",this_buildcomp,".limma_voom_dream.all.txt.gz")
      mydf_anylogfc<-fread(my_filename_anylogfc)%>%mutate(tissue=this_tissue, comparison=this_buildcomp)
      all_diffexp_anylogfc=rbind(all_diffexp_anylogfc,mydf_anylogfc)
      
    }
    
  }
  all_diffexp_anylogfc[, c("t","P.Value","z.std","AveExpr"):=NULL] 
  all_diffexp_anylogfc_casted=dcast(all_diffexp_anylogfc,gene+tissue~comparison,value.var=c('logFC','adj.P.Val'))
  return(all_diffexp_anylogfc_casted)
}
get_ave_zchange<-function(){}
get_chm13_dic<-function(){
  convert_chm13=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/chm13_enstANDensg.txt")
  convert_chm13[,c("V2","V4","V5"):=NULL]
  colnames(convert_chm13)=c("ensg","chm13_gene")
  convert_chm13=unique(convert_chm13)
  chm13_gene_dic=convert_chm13%>%mutate(ensg=gsub("\\..*","",ensg))%>%pull(ensg)
  names(chm13_gene_dic)<-convert_chm13$chm13_gene
  return(chm13_gene_dic)
}
get_mapping<-function(chm13_gene_dic){
  mytissues=c("Blood","Fibroblast","Muscle","PBMC","iPSCNPC","iPS")
  all_mapping=data.frame()
  for(this_tissue in mytissues){
    print(this_tissue)
    for(this_build in c("hg19","hg38","chm13")){
      print(this_build)
      my_filename=paste0("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/Mapping/mapping_",this_build,"_",this_tissue,"_allsamples_totals_meangene.txt.gz")
      if(this_build !="chm13"){
        mydf<-fread(my_filename)%>%mutate(tissue=this_tissue, build=this_build)%>%
          mutate(mean_prop_multimapped=mean_multi_mapped/(mean_multi_mapped+mean_uniquely_mapped))%>%
          mutate(ensg=gsub("\\..*","",ensg))
      }else{

        mydf<-fread(my_filename)%>%mutate(tissue=this_tissue, build=this_build)%>%
          mutate(ensg=chm13_gene_dic[chm13_gene])%>%
          mutate(mean_prop_multimapped=mean_multi_mapped/(mean_multi_mapped+mean_uniquely_mapped))
        mydf=mydf[,'chm13_gene':=NULL]
      }
     
      all_mapping=rbind(all_mapping,mydf)
    }
  }
  all_mapping=all_mapping%>%dplyr::filter(!is.na(ensg)&ensg!="")%>%unique()
  ###NOTE: IF THERE ARE MULTIPLE CHM13 TRANSCRIPTS DETECTED, THE MEAN IS TAKEN OF THEM FOR NUM TRANSCRIPTS AND PROP MULTIMAPPED
  all_mapping_casted=dcast(all_mapping,ensg+tissue~build,value.var=c('mean_prop_multimapped','mean_num_transcripts'),fun.aggregate=mean)
  
  all_mapping_casted$diff_mean_propmultimapped_hg19_hg38=all_mapping_casted$mean_prop_multimapped_hg38-all_mapping_casted$mean_prop_multimapped_hg19
  all_mapping_casted$diff_mean_propmultimapped_hg38_chm13=all_mapping_casted$mean_prop_multimapped_chm13-all_mapping_casted$mean_prop_multimapped_hg38
  all_mapping_casted$diff_num_transcripts_hg19_hg38=all_mapping_casted$mean_num_transcripts_hg38-all_mapping_casted$mean_num_transcripts_hg19
  all_mapping_casted$diff_num_transcripts_hg38_chm13=all_mapping_casted$mean_num_transcripts_chm13-all_mapping_casted$mean_num_transcripts_hg38
  return(all_mapping_casted)
}
get_omim<-function(min_confidence){
  library(biomaRt)
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
                           ifelse(confidence == 2, "linkage_mapping", 
                                  ifelse(confidence == 3, "known_molec_cause", 
                                         ifelse(confidence == 4, "del_dup_syndrome", NA))))) %>% 
    separate_rows(gene, convert = T, sep=", ") %>% 
    mutate(ensg = mim_to_ensg[mim]) %>% 
    mutate(ensg=ifelse(ensg=="" | is.na(ensg),hgnc_to_ensg[gene],ensg))%>%
    mutate(ensg=ifelse(ensg=="" | is.na(ensg),ensg_to_hgnc_more[gene],ensg))%>%
    
    distinct() 
  pheno_filtered<-pheno_clean%>%dplyr::filter(confidence>=min_confidence)
  return(pheno_filtered)
}
get_annotations<-function(){
  cosmic=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/cosmic_cancer_genes.tsv")
  cosmic$ensg=gsub("\\..*","",str_extract(cosmic$Synonyms,"ENSG(.*)\\."))
  cosmic=cosmic[,c("Somatic","Germline","Cancer Syndrome","ensg")]%>%dplyr::filter(!is.na(ensg))
  
  opentargets=fread("/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/opentargetscores.txt")
  colnames(opentargets)<-c("target_id","ensg","target_score","x")
  opentargets_summ=opentargets%>%group_by(ensg)%>% dplyr::filter(target_score ==max(target_score)) %>%dplyr::select(-x) #
  
  cosmic_opentargets=merge(cosmic,opentargets_summ,by="ensg",all=T)
  
  omim_df=get_omim(min_confidence=0)
  omim=omim_df%>%dplyr::select(ensg,confidence,pheno)%>%dplyr::filter(!is.na(ensg))%>%dplyr::filter(ensg!='')%>%unique()
  cosmic_opentargets_omim=merge(cosmic_opentargets,omim,all=T)
  
  hgnc_all=fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/hgnc_ensg_convert_complete_annotated.txt")
  hgnc_all=hgnc_all%>%dplyr::select(-strand,-n_TSS,-gene_start,-gene_end)
  cosmic_opentargets_omim_hgnc=merge(cosmic_opentargets_omim,hgnc_all,by="ensg",all=TRUE)
  return(cosmic_opentargets_omim_hgnc%>%mutate(omim_confidence=confidence)%>%dplyr::select(-confidence))
}

main<-function(){
  diffexp=get_diffexp()
  chm13_gene_dic=get_chm13_dic()
  mapping_stats=get_mapping(chm13_gene_dic)
  diffexp_mapping=merge(diffexp,mapping_stats,by.x=c('gene','tissue'),by.y=c('ensg','tissue'),all.x = TRUE)%>%mutate(ensg=gene)%>%dplyr::select(-gene)
  geneannotations=get_annotations()
  diffexp_mapping_annotated=merge(diffexp_mapping,geneannotations,all.x=T,by="ensg")%>%dplyr::filter(ensg!="None")
  fwrite(diffexp_mapping_annotated,file="/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/final_diffquant_table.txt.gz",quote=FALSE,sep="\t")
}

main()