library(data.table)
library(ggpubr)
library("biomaRt")
library('tidyverse')

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
opentargets_summ=opentargets%>%group_by(ensg)%>% dplyr::filter(target_score ==max(target_score)) %>%select(-x) #

#TPM
read_in_file<-function(filename){
  my_df=fread(filename)
  colnames(my_df)=c("sample","gene","transcript","tpm")
  my_df$gene<-sapply(strsplit(my_df$gene,"\\."),"[[",1)
  return(my_df)
}
get_build_specific_expression<-function(build1_tissue,build2_tissue,notexp,build1_name,build2_name,this_tiss){
  notexp_summed=notexp%>%mutate(single_build_exp=(as.numeric(get(paste0("exp_",build1_name)))+
                                                     as.numeric(get(paste0("exp_",build2_name)))))
  exp_one_build=notexp_summed%>%dplyr::filter(single_build_exp==1)%>%
    mutate(which_build_is_expressed=ifelse(get(paste0("exp_",build1_name))==TRUE,build1_name,build2_name))
  
  build1_tissue_filt=build1_tissue%>%dplyr::filter(gene %in%exp_one_build$gene)
  build2_tissue_filt=build2_tissue%>%dplyr::filter(gene %in%exp_one_build$gene)
  
  build1_summary=build1_tissue_filt[,as.list(quantile(tpm, c(.25, .5, .75)), na.rm = TRUE), by = gene]
  colnames(build1_summary)[2:4]= c(paste0("Q1_build1"),paste0("Q2_build1"),paste0("Q3_build1"))
  build2_summary=build2_tissue_filt[,as.list(quantile(tpm, c(.25, .5, .75)), na.rm = TRUE), by = gene]
  colnames(build2_summary)[2:4]= c(paste0("Q1_build2"),paste0("Q2_build2"),paste0("Q3_build2"))
  summaries_merged<-merge(build1_summary,build2_summary,by="gene")
  summaries_merged$expressed_build_med= pmax(summaries_merged[,3],summaries_merged[,6])
  all_merged<-merge(exp_one_build[,c("gene","which_build_is_expressed")],summaries_merged,by="gene")
  
  diff_exp_return=all_merged%>%
    mutate(comparisons=paste0(build1_name,"-",build2_name))%>%
    mutate(tissue=this_tiss)
  return(diff_exp_return)
}
all_build_specific_expression<-function(){
  ##files to read in
  #expressed or not
  chm13_hg38_blood_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/Blood.chm13ensembl_vs_hg38.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")%>%dplyr::filter(!(gene %in% chrY_chm13_genes$gene))
  chm13_hg38_fibroblast_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/Fibroblast.chm13ensembl_vs_hg38.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")%>%dplyr::filter(!(gene %in% chrY_chm13_genes$gene))
  chm13_hg38_iPS_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/iPS.chm13ensembl_vs_hg38.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")%>%dplyr::filter(!(gene %in% chrY_chm13_genes$gene))
  chm13_hg38_iPSC_NPC_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/iPSC_NPC.chm13ensembl_vs_hg38.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")%>%dplyr::filter(!(gene %in% chrY_chm13_genes$gene))
  chm13_hg38_muscle_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/Muscle.chm13ensembl_vs_hg38.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")%>%dplyr::filter(!(gene %in% chrY_chm13_genes$gene))
  chm13_hg38_PBMC_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/PBMC.chm13ensembl_vs_hg38.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")%>%dplyr::filter(!(gene %in% chrY_chm13_genes$gene))
  
  hg38_hg19_blood_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/Blood.hg38_vs_hg19.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")
  hg38_hg19_fibroblast_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/Fibroblast.hg38_vs_hg19.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")
  hg38_hg19_iPS_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/iPS.hg38_vs_hg19.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")
  hg38_hg19_iPSC_NPC_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/iPSC_NPC.hg38_vs_hg19.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")
  hg38_hg19_muscle_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/Muscle.hg38_vs_hg19.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")
  hg38_hg19_PBMC_notexp=fread("/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/02_intermediateOutputs/PBMC.hg38_vs_hg19.dedupOptical_minMQ255_rsem.genes.not_enough_exp.0counts_0.1cpm_30pct_samples.txt")
  
  
  hg19_fibroblast=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Fibroblast.hg19.dedupOptical_minMQ255_rsem.genes.results")
  hg19_blood=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Blood.hg19.dedupOptical_minMQ255_rsem.genes.results")
  hg19_iPS=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/iPS.hg19.dedupOptical_minMQ255_rsem.genes.results")
  hg19_iPSC_NPC=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/iPSC_NPC.hg19.dedupOptical_minMQ255_rsem.genes.results")
  hg19_muscle=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/Muscle.hg19.dedupOptical_minMQ255_rsem.genes.results")
  hg19_PBMC=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg19/eOutliers/PBMC.hg19.dedupOptical_minMQ255_rsem.genes.results")
  
  hg38_fibroblast=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Fibroblast.hg38.dedupOptical_minMQ255_rsem.genes.results")
  hg38_blood=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Blood.hg38.dedupOptical_minMQ255_rsem.genes.results")
  hg38_iPS=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/iPS.hg38.dedupOptical_minMQ255_rsem.genes.results")
  hg38_iPSC_NPC=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/iPSC_NPC.hg38.dedupOptical_minMQ255_rsem.genes.results")
  hg38_muscle=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/Muscle.hg38.dedupOptical_minMQ255_rsem.genes.results")
  hg38_PBMC=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/PBMC.hg38.dedupOptical_minMQ255_rsem.genes.results")
  
  chm13_fibroblast=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Fibroblast.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results")
  chm13_blood=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Blood.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results")
  chm13_iPS=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/iPS.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results")
  chm13_iPSC_NPC=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/iPSC_NPC.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results")
  chm13_muscle=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/Muscle.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results")
  chm13_PBMC=read_in_file("/oak/stanford/groups/smontgom/shared/UDN/Output/chm13/eOutliers/PBMC.chm13ensembl.dedupOptical_minMQ255_rsem.genes.results")
  
  hg19_hg38_fibroblast<-get_build_specific_expression(hg19_fibroblast,hg38_fibroblast,hg38_hg19_fibroblast_notexp,"hg19","hg38","fibroblast")
  hg19_hg38_blood<-get_build_specific_expression(hg19_blood,hg38_blood,hg38_hg19_blood_notexp,"hg19","hg38","blood")
  hg19_hg38_iPS<-get_build_specific_expression(hg19_iPS,hg38_iPS,hg38_hg19_iPS_notexp,"hg19","hg38","iPS")
  hg19_hg38_iPSC_NPC<-get_build_specific_expression(hg19_iPSC_NPC,hg38_iPSC_NPC,hg38_hg19_iPSC_NPC_notexp,"hg19","hg38","iPSC_NPC")
  hg19_hg38_muscle<-get_build_specific_expression(hg19_muscle,hg38_muscle,hg38_hg19_muscle_notexp,"hg19","hg38","muscle")
  hg19_hg38_PBMC<-get_build_specific_expression(hg19_PBMC,hg38_PBMC,hg38_hg19_PBMC_notexp,"hg19","hg38","PBMC")
  
  hg38_chm13_fibroblast<-get_build_specific_expression(hg38_fibroblast,chm13_fibroblast,chm13_hg38_fibroblast_notexp,"hg38","chm13ensembl","fibroblast")
  hg38_chm13_blood<-get_build_specific_expression(hg38_blood,chm13_blood,chm13_hg38_blood_notexp,"hg38","chm13ensembl","blood")
  hg38_chm13_iPS<-get_build_specific_expression(hg38_iPS,chm13_iPS,chm13_hg38_iPS_notexp,"hg38","chm13ensembl","iPS")
  hg38_chm13_iPSC_NPC<-get_build_specific_expression(hg38_iPSC_NPC,chm13_iPSC_NPC,chm13_hg38_iPSC_NPC_notexp,"hg38","chm13ensembl","iPSC_NPC")
  hg38_chm13_muscle<-get_build_specific_expression(hg38_muscle,chm13_muscle,chm13_hg38_muscle_notexp,"hg38","chm13ensembl","muscle")
  hg38_chm13_PBMC<-get_build_specific_expression(hg38_PBMC,chm13_PBMC,chm13_hg38_PBMC_notexp,"hg38","chm13ensembl","PBMC")
  
  all_comparisons=rbind(hg19_hg38_fibroblast,
                        hg19_hg38_blood,
                        hg19_hg38_iPS,
                        hg19_hg38_iPSC_NPC,
                        hg19_hg38_muscle,
                        hg19_hg38_PBMC,
                        hg38_chm13_fibroblast,
                        hg38_chm13_blood,
                        hg38_chm13_iPS,
                        hg38_chm13_iPSC_NPC,
                        hg38_chm13_muscle,
                        hg38_chm13_PBMC) 
  return(all_comparisons)
}
all_comparisons=all_build_specific_expression()
#fwrite(all_comparisons,"/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/all_build_specific_expression.txt.gz")

BMS1P8_blood<-rbind(cbind(hg19_blood%>%dplyr::filter(gene=="ENSG00000260518"),build='hg19'),
              cbind(hg38_blood%>%dplyr::filter(gene=="ENSG00000260518"),build='hg38'),
              cbind(chm13_blood%>%dplyr::filter(gene=="ENSG00000260518"),build='chm13'))
BMS1P8_blood$build=factor(BMS1P8_blood$build,levels=c('hg19','hg38','chm13'))
BMS1P8_fibroblast<-rbind(cbind(hg19_fibroblast%>%dplyr::filter(gene=="ENSG00000260518"),build='hg19'),
                    cbind(hg38_fibroblast%>%dplyr::filter(gene=="ENSG00000260518"),build='hg38'),
                    cbind(chm13_fibroblast%>%dplyr::filter(gene=="ENSG00000260518"),build='chm13'))
BMS1P8_fibroblast$build=factor(BMS1P8_fibroblast$build,levels=c('hg19','hg38','chm13'))


ggplot(BMS1P8_fibroblast,aes(x=build,y=tpm,fill=build))+
  geom_boxplot()+
  geom_jitter(width = 0.15)+
  ggtitle('Fibroblast BMS1P8')+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  scale_fill_manual(values=c( "#8C1515","#dea700","#01427a"))+
  theme_bw(base_size = 20)

opentargets_summ0.6=opentargets_summ%>%dplyr::filter(target_score>0.6); opentargets_summ0.8=opentargets_summ%>%dplyr::filter(target_score>0.8)
table(unique(all_comparisons$gene)%in%opentargets_summ0.6$ensg); length(unique(opentargets_summ0.6$ensg));table(unique(all_comparisons$gene)%in%opentargets_summ0.8$ensg); length(unique(opentargets_summ0.8$ensg))
(table(unique(all_comparisons$gene)%in%unique(omim_genes))); print(length(unique(omim_genes)))
comp_opentarget0.8=all_comparisons%>%dplyr::filter(gene %in% opentargets_summ0.8$ensg)
opentargets%>%dplyr::filter(ensg %in% comp_opentarget0.8$gene)%>%arrange(desc(target_score))

all_comparisons_anno=merge(all_comparisons,opentargets_summ0.8,by.x = 'gene',by.y = "ensg",all.x = T)
all_comparisons_anno=merge(all_comparisons_anno,omim_df[,c('ensg','gene','confidence')],by.x = 'gene',by.y = "ensg",all.x = T)%>%unique()

distr_plot=ggplot(all_comparisons,aes(x=expressed_build_med,color=tissue,alpha=0.5))+
  geom_density()+
  scale_x_continuous(trans="log10",labels=scales::comma)+
  theme_bw()+
  scale_color_manual(values=c("#D90025","#27A795","#0470B5","#BD71DC","#AACA2F","#EF7C18"))+
  xlab("TPM expressed build")+
  facet_wrap(~comparisons)

number_genes_exp_1build=as.data.frame(rbind(c(comparison="hg19-hg38",tissue="fibroblast",num_build_specific_exp=nrow(hg19_hg38_fibroblast)),
                                            c(comparison="hg19-hg38",tissue="blood",num_build_specific_exp=nrow(hg19_hg38_blood)),
                                            c(comparison="hg19-hg38",tissue="iPS",num_build_specific_exp=nrow(hg19_hg38_iPS)),
                                            c(comparison="hg19-hg38",tissue="iPSC_NPC",num_build_specific_exp=nrow(hg19_hg38_iPSC_NPC)),
                                            c(comparison="hg19-hg38",tissue="muscle",num_build_specific_exp=nrow(hg19_hg38_muscle)),
                                            c(comparison="hg19-hg38",tissue="PBMC",num_build_specific_exp=nrow(hg19_hg38_PBMC)),
                                            c(comparison="hg38-chm13",tissue="fibroblast",num_build_specific_exp=nrow(hg38_chm13_fibroblast)),
                                            c(comparison="hg38-chm13",tissue="blood",num_build_specific_exp=nrow(hg38_chm13_blood)),
                                            c(comparison="hg38-chm13",tissue="iPS",num_build_specific_exp=nrow(hg38_chm13_iPS)),
                                            c(comparison="hg38-chm13",tissue="iPSC_NPC",num_build_specific_exp=nrow(hg38_chm13_iPSC_NPC)),
                                            c(comparison="hg38-chm13",tissue="muscle",num_build_specific_exp=nrow(hg38_chm13_muscle)),
                                            c(comparison="hg38-chm13",tissue="PBMC",num_build_specific_exp=nrow(hg38_chm13_PBMC))))
number_genes_exp_1build$num_build_specific_exp<-as.numeric(as.character(number_genes_exp_1build$num_build_specific_exp))
barplot=ggplot(number_genes_exp_1build,aes(y=num_build_specific_exp,x=tissue,fill=tissue))+
  geom_bar(stat="identity",position = "dodge")+coord_flip()+
  facet_wrap(~comparison)+
  scale_fill_manual(values=c("#D90025","#27A795","#0470B5","#BD71DC","#AACA2F","#EF7C18"))+
  theme_bw()

ggarrange(distr_plot,barplot, 
  nrow = 2)

hg19hg38_df=all_comparisons%>%dplyr::filter(comparisons=="hg19-hg38")%>%group_by(tissue)
hg38chm13ensembl_df=all_comparisons%>%dplyr::filter(comparisons=="hg38-chm13ensembl")%>%group_by(tissue)
in_one_not_other=merge(hg19hg38_df%>%summarise(in_hg19_not_hg38=sum(which_build_is_expressed=="hg19"),
                        not_hg19_in_hg38=sum(which_build_is_expressed=="hg38")),
                  hg38chm13ensembl_df%>%summarise(in_hg38_not_chm13ensembl=sum(which_build_is_expressed=="hg38"),
                                                  not_hg38_in_chm13ensembl=sum(which_build_is_expressed=="chm13ensembl")))
median_medians=as.data.table(all_comparisons)[,(as.list(quantile(expressed_build_med, c(0,.25, .5, .75,1)), na.rm = TRUE)), 
                                              by = c("tissue","which_build_is_expressed","comparisons")]
median_medians_genes=as.data.table(all_comparisons)[,.N, by = c("tissue","which_build_is_expressed","comparisons")]



exp_with_omim=all_tpms%>%dplyr::filter(omim=="omim")%>%pull(gene)%>%unique()
opentargets_summ0.6=opentargets_summ%>%dplyr::filter(target_score>0.6); opentargets_summ0.8=opentargets_summ%>%dplyr::filter(target_score>0.8)
table(unique(all_diffexp$gene)%in%opentargets_summ0.6$ensg); length(unique(opentargets_summ0.6$ensg));table(unique(all_diffexp$gene)%in%opentargets_summ0.8$ensg); length(unique(opentargets_summ0.8$ensg))
length(table(unique(all_diffexp$gene)%in%unique(omim_genes))); print(length(unique(omim_genes)))

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
ggplot(hg38chm13ensembl_df,aes(x=tissue,y=expressed_build_med,fill=tissue,alpha=which_build_is_expressed))+
  geom_split_violin(aes(group=interaction(which_build_is_expressed,tissue)))+
  scale_fill_manual(values=c("#D90025","#27A795","#0470B5","#BD71DC","#AACA2F","#EF7C18"))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  ggtitle("expressed in just one build")+
  scale_alpha_discrete(range = c(0.3, 0.9), guide = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw()
