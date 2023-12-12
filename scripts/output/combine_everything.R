library(data.table)
library(optparse)
library(dplyr)

#--- OPTION PARSER
option_list = list(
  make_option(c( "--tpm"), type="character", default=NULL, help="tpm file", metavar="character"),
  make_option(c( "--sample"), type="character", default=NULL, help="tpm file", metavar="character"),
  make_option(c("--med_counts"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--vcf"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--build"), type="character", default=NULL, help="build", metavar="character"),
  make_option(c("--genes_or_isoforms"), type="character", default=NULL, help="genes_or_isoforms", metavar="character"),
  make_option(c("--expression_zscore_udn"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--expression_zscore_wgtex"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--splicing_zscore"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--ase_zscore"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--anno_file"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--metadata"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--outfile"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--outfile_outliers_spl"), type="character", default=NULL, help="median tpm file", metavar="character"),
  make_option(c("--outfile_outliers_exp"), type="character", default=NULL, help="median tpm file", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
tpm_file<-as.character(opt$tpm)
med_counts_file<-as.character(opt$med_counts)
build<-as.character(opt$build)
genes_or_isoforms<-as.character(opt$genes_or_isoforms)
vcf_file<-as.character(opt$vcf)
anno_file<-as.character(opt$anno_file)
expression_zscore_udn_file<-as.character(opt$expression_zscore_udn)
expression_zscore_wgtex_file<-as.character(opt$expression_zscore_wgtex)
splicing_zscore_file<-as.character(opt$splicing_zscore)
ase_zscore_file<-as.character(opt$ase_zscore)
metadata_file<-as.character(opt$metadata)
sample<-as.character(opt$sample)
outfile<-as.character(opt$outfile)
outfile_exp<-as.character(opt$outfile_outliers_exp)
outfile_spl<-as.character(opt$outfile_outliers_spl)
# 
#  sample="XXX"
# tpm_file<-"/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/RSEM/XXX.hg38.sorted_opticalDup100_dedupOptical_minMQ255_rsem.genes.results"
#  med_counts_file<-"/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparison/TPMDifferences/Fibroblast_tpm_genesdiff.txt.gz"
# vcf_file<-"/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/VCFs/RVAnnotated/YYY_cadd.tsv.gz"
# expression_zscore_udn_file<-"/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/eOutliers/SplitByInd/XXX.genes.txt"
# # expression_zscore_wgtex_file<-"NoFile"
# splicing_zscore_file<-"/oak/stanford/groups/smontgom/shared/UDN/Output/hg38/Splicing/SplitByInd/XXX_genes.hg38.sorted_opticalDup100_dedupOptical_minMQ255.collapsed.bed"
# ase_zscore_file<-"NoFile"
# build<-"hg38"
# genes_or_isoforms<-"genes"
# metadata_file<-"/oak/stanford/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv"
# anno_file<-"/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/combined_gene_disease_info.txt.gz"
# outfile<-"/oak/stanford/groups/smontgom/shared/UDN/PreprocessingHG38Primary/Combined/XXX.txt"

md<-fread(metadata_file)
md_sample<-md%>%select(-seq_batch)%>% dplyr::filter(sample_id==sample) #remove duplicate name seq_batch need to fix md
udn_id<-md_sample$indv_id
sample_hpo<-gsub(",",";",md_sample$HPO_terms_ID)


#read in, subset, and rename gene level fiels
gene_anno<-fread(anno_file)
tpm<-fread(tpm_file)
print("TPM")
if(genes_or_isoforms=="genes"){
  colnames(tpm)[1:2]<-c("ensg","enst")
  leveltype="ensg"
  # colnames(med_counts_all)[c(1)]<-c("ensg")
  
}else if(genes_or_isoforms=="isoforms"){
  colnames(tpm)[1:2]<-c("enst","ensg")
  # colnames(med_counts_all)[c(1)]<-c("enst")
  leveltype="enst"
  
}else{
  stop("GENES_OR_ISOFORMS IS NOT DEFINED")
}
tpm$ensg<-sapply(strsplit(tpm$ensg,"\\."),"[[",1)
tpm$enst<-sapply(strsplit(tpm$enst,"\\."),"[[",1)

#join gene level information with TPMs
tpm_red<-tpm[,c("ensg","enst","TPM")]

if(med_counts_file=="NoFile"){
  tpm_gtexmed<-cbind(tpm_red,med_counts="NA")
}else{
  med_counts_all<-fread(med_counts_file)
  med_counts<-data.frame("id"=sapply(strsplit(unlist(med_counts_all[,1]),"\\."),"[[",1),
                         "med_counts"= med_counts_all %>% pull(paste0("med_",build)))
  tpm_gtexmed<-merge(tpm_red,med_counts,by.x=leveltype,by.y="id",all.x=T)
  
}



# med_counts$ensg<-sapply(strsplit(med_counts$ensg,"\\."),"[[",1)



#join gtex med counts with tpms
tpm_gtexmed_anno<-merge(tpm_gtexmed,gene_anno,by="ensg",all.x=T)
#sample_hpo
print("exp")
##sample_hpo_vec<-unlist(strsplit(sample_hpo,";"))
#tpm_gtexmed_anno$clinvar_hpo_num_intersects<-unlist(lapply(tpm_gtexmed_anno %>% pull("clinvar_hpo"),function(x){
#       length(intersect(unlist(strsplit(x,";")),sample_hpo_vec))}))
#tpm_gtexmed_anno$DDG_hpo_num_intersects<-unlist(lapply(tpm_gtexmed_anno %>% pull("DDG_hpo_terms"),function(x){
#  length(intersect(unlist(strsplit(x,";")),sample_hpo_vec))}))


##check if expression information, if so add, if not DON'T
if(expression_zscore_udn_file =="NoFile"){
  tpm_gtexmed_anno_exp<-cbind(tpm_gtexmed_anno,"z_exp_udn"="NA","rank_exp_all"="NA","rank_exp_under"="NA","rank_exp_over"="NA")
  # tpm_gtexmed_anno_exp<-data.frame("z_exp_udn"="NA","rank_exp_all"="NA","rank_exp_under"="NA","rank_exp_over"="NA")
}else{
  expression_zscore_udn<-fread(expression_zscore_udn_file)
  colnames(expression_zscore_udn)<-c("sample_id",leveltype,"z_exp_udn","rank_exp_all","rank_exp_under","rank_exp_over")
  expression_zscore_udn=expression_zscore_udn%>%mutate(rank_exp_all=rank(-abs(z_exp_udn),na.last=TRUE,))  %>%
                         mutate(rank_exp_over=rank(-(z_exp_udn),na.last=TRUE,))%>%
                         mutate(rank_exp_under=rank((z_exp_udn),na.last=TRUE,))
  expression_zscore_udn<-expression_zscore_udn %>%mutate(!!leveltype:=sapply(strsplit(expression_zscore_udn %>% pull(leveltype),"\\."),"[[",1))
  tpm_gtexmed_anno_exp<-merge(tpm_gtexmed_anno,expression_zscore_udn,by=leveltype,all.x=T)
   # tpm_gtexmed_anno_exp<-expression_zscore_udn[,2:3]
  # tpm_gtexmed_anno_exp<-tpm_gtexmed_anno[expression_zscore_udn[,2:3],on="ensg"]
  
}
print("exp w gtex")
#add z score exp w gtex
if(expression_zscore_wgtex_file =="NoFile"){
  print("NO FILE FOUNDS GTEX ")
  tpm_gtexmed_anno_exp2<-tpm_gtexmed_anno_exp
  tpm_gtexmed_anno_exp2$z_exp_udnWITHgtex<-"NA"
  
}else{
  expression_zscore_wgtex<-fread(expression_zscore_wgtex_file)
  colnames(expression_zscore_wgtex)<-c("sample_id","ensg","z_exp_udnWITHgtex")
  print(head(expression_zscore_wgtex))
  print("ensg")
  expression_zscore_wgtex$ensg<-sapply(strsplit(expression_zscore_wgtex$ensg,"\\."),"[[",1)
  print(head(expression_zscore_wgtex[,2:3]))
  tpm_gtexmed_anno_exp2<-merge(tpm_gtexmed_anno_exp,expression_zscore_wgtex[,2:3],by="ensg",all=T)
  # tpm_gtexmed_anno_exp2<-tpm_gtexmed_anno_exp[expression_zscore_wgtex[,2:3],on="ensg"]
}
print("what happened above?")
print(head(tpm_gtexmed_anno_exp2))
print("spl")

##check if splicing information, if so add, if not DON'T
if(splicing_zscore_file =="NoFile"){
  tpm_gtexmed_anno_exp2_spl<-tpm_gtexmed_anno_exp2
  tpm_gtexmed_anno_exp2_spl$chr<-NA
  tpm_gtexmed_anno_exp2_spl$exon_junc_start<-NA
  tpm_gtexmed_anno_exp2_spl$exon_junc_end<-NA
  tpm_gtexmed_anno_exp2_spl$cluster_id<-NA
  tpm_gtexmed_anno_exp2_spl$sample_id<-NA
  tpm_gtexmed_anno_exp2_spl$z_spl<-NA 
  tpm_gtexmed_anno_exp2_spl$rank_overall<-NA
  tpm_gtexmed_anno_exp2_spl$rank_withinind<-NA
  tpm_gtexmed_anno_exp2_spl$rank_withincluster<-NA
  tpm_gtexmed_anno_exp2_spl$junction_id<-NA
  tpm_gtexmed_anno_exp2_spl$num_reads_supporting_junc<-NA
  tpm_gtexmed_anno_exp2_spl$strand<-NA
  tpm_gtexmed_anno_exp2_spl$splice_site<-NA
  tpm_gtexmed_anno_exp2_spl$acceptors_skipped<-NA
  tpm_gtexmed_anno_exp2_spl$exons_skipped<-NA
  tpm_gtexmed_anno_exp2_spl$donors_skipped<-NA  
  tpm_gtexmed_anno_exp2_spl$anchor<-NA
  tpm_gtexmed_anno_exp2_spl$known_donor<-NA
  tpm_gtexmed_anno_exp2_spl$known_junction<-NA
  # tpm_gtexmed_anno_exp2_spl$enst_spl<-NA
  # tpm_gtexmed_anno_exp2_spl$ensg<-NA
}else{
  splicing_zscore<-fread(splicing_zscore_file)
  print(head(splicing_zscore))
  colnames(splicing_zscore)<-c("chr","exon_junc_start","exon_junc_end","cluster_id","sample_id",
                               "z_spl","junction_id","num_reads_supporting_junc",
                               "strand","splice_site","acceptors_skipped","exons_skipped",
                               "donors_skipped","anchor","known_donor","known_junction","enst",
                               "ensg","na") #,"enst")
  print(head(splicing_zscore))
  splicing_zscore$ensg<-sapply(strsplit(splicing_zscore$ensg,"\\."),"[[",1)
  splicing_zscore=splicing_zscore%>%mutate(rank_withinind=rank(-(z_spl),na.last=TRUE,))
  #for each group (gene), reduce by choosing max z score!
  print("Now head redbyensg")
  splicing_zscore_redbyensg<-setDT((splicing_zscore))[, .SD[which.max((z_spl))], by=leveltype]
  print(head(splicing_zscore_redbyensg))
  tpm_gtexmed_anno_exp2_spl<-merge(tpm_gtexmed_anno_exp2,splicing_zscore_redbyensg[,-"na"],by=leveltype,all.x=T,allow.cartesian = T)
}
print(head(tpm_gtexmed_anno_exp2_spl))
print("ase")

###ase
if(ase_zscore_file =="NoFile"){
  tpm_gtexmed_anno_exp2_spl_ase<-tpm_gtexmed_anno_exp2_spl
  tpm_gtexmed_anno_exp2_spl_ase$z_ase<-"NA"
}else{
  #FILL IN 
}
print(head(tpm_gtexmed_anno_exp2_spl_ase))
print("vcf")

###vcf
if(vcf_file=="NoFile"){
  print("NO VCF")
  tpm_gtexmed_anno_exp2_spl_ase_vcf<-tpm_gtexmed_anno_exp2_spl_ase
  tpm_gtexmed_anno_exp2_spl_ase_vcf$chr<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$pos<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$ref<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$alt<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$ind<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$gene<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$enst_variant<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$maf<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$cadd_phred<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$impactful<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$impact_score<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$vartype<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$dna_change<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$aa_change<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$disease<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$gene_fam<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$spl_score<-NA
  tpm_gtexmed_anno_exp2_spl_ase_vcf$num_rv<-NA

}else{
  vcf<-fread(vcf_file)
  colnames(vcf)<-c(colnames(vcf)[1:13],"dna_change","aa_change","disease","gene_fam","spl_score","num_rv")
  colnames(vcf)[1:4]<-paste0("rv_",colnames(vcf)[1:4])
  colnames(vcf)[8]<-"enst"
  # vcf_subset<-vcf[,c("chr","pos","ref","alt","gene","ensg","maf","cadd_phred","impactful","impact_score",
  #                    "vartype","num_rv","dna_change","aa_change")]
  # tpm_gtexmed_anno_exp2_spl_ase_vcf<-tpm_gtexmed_anno_exp2_spl_ase[vcf,on="ensg"]
  print("HELLLOOO")
  print(head(tpm_gtexmed_anno_exp2_spl_ase))
  print("now vcf")
  print(head(vcf))
  tpm_gtexmed_anno_exp2_spl_ase_vcf<-merge(tpm_gtexmed_anno_exp2_spl_ase,vcf,by=leveltype,all.x=T)
  print("AFTER ????")
}
print("sort")

# sorted_output<-tpm_gtexmed_anno_exp2_spl_ase_vcf%>%arrange(impact_score,cadd_phred,z_exp_udn,z_spl)
#sorted_output_simplified<-tpm_gtexmed_anno_exp2_spl_ase_vcf %>% select(-hgnc_id,-ncbi_id,-gene_fam,-ind,-gene,-disease,-spl_score) %>%
sorted_output_simplified<-tpm_gtexmed_anno_exp2_spl_ase_vcf %>% select(-hgnc_id,-ncbi_id,-gene_fam,-ind,-gene,-disease,-spl_score) %>%
  rename(cytoband=chromosome,z_splice=z_spl,#rank_splice=rank_spl,#splice_transcript=enst_spl,
         gene_family=gene_group,
	 DDG_disease=DDG,clinvar_disease=clinvar,
         #rv_chr=chr, rv_pos=pos,rv_ref=ref,rv_alt=alt, #rv_transcript=enst_variant,
         rv_maf=maf,rv_cadd_phred=cadd_phred,is_impactful=impactful, rv_impact_score=impact_score,rv_vartype=vartype,
         rv_dna_change=dna_change,rv_aa_change=aa_change,rv_num=num_rv,
	 splice_rank_withinind=rank_withinind)%>%
  arrange(rv_impact_score,rv_cadd_phred,z_exp_udn,z_splice)
fwrite(sorted_output_simplified,file=outfile,sep="\t",quote=F,na = NA)

print("num rv")
print(nrow(sorted_output_simplified%>%dplyr::filter(!is.na(rv_num))))

sorted_output_simplified_sOutliersRVsONLY=sorted_output_simplified%>%
  dplyr::filter(!is.na(rv_num)&(abs(z_splice)>2))
fwrite(sorted_output_simplified_sOutliersRVsONLY,file=outfile_spl,sep="\t",quote=F,na = NA)

sorted_output_simplified_eOutliersRVsONLY=sorted_output_simplified%>%
  dplyr::filter(!is.na(rv_num)&(abs(z_exp_udn)>2))
fwrite(sorted_output_simplified_eOutliersRVsONLY,file=outfile_exp,sep="\t",quote=F,na = NA)



