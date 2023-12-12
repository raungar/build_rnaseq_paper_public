library(data.table)
library(dplyr)
library(optparse)

# infile="/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/Mapping/mapping_hg19_allsamples_summarised.txt.gz"
# gtf_file="/Volumes/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/ReferenceFiles/hg19_enstANDensg.txt"
# build="hg19"

parse_my_args<-function(){
  option_list = list(
    make_option(c( "--infile"), type="character", default=NULL, help="z hg19", metavar="character"),
    make_option(c( "--gtf_file"), type="character", default=NULL, help="z hg19", metavar="character"),
    make_option(c( "--build"), type="character", default=NULL, help="z hg19", metavar="character"),
    make_option(c("--outfile"), type="character", default=NULL, help="outfile", metavar="character"),
    make_option(c("--outfile_mean"), type="character", default=NULL, help="outfile", metavar="character")
  )
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt)
}
get_gtf<-function(opt){
  gtf=fread(opt$gtf_file,header=F)
  setnames(gtf,c("ensg","enst"))
  gtf_dic=gtf$ensg
  names(gtf_dic)<-gtf$enst
  return(gtf_dic)
}
get_gtf_chm13<-function(opt){
  gtf=fread(opt$gtf_file,header=F)
  setnames(gtf,c("ensg","enst","chm13_gene","chm13_transcript","na"))
  gtf_chm13_to_ensembl_dic=gtf$ensg
  names(gtf_chm13_to_ensembl_dic)=gtf$chm13_gene
  gtf_dic=gtf$chm13_gene
  names(gtf_dic)=gtf$chm13_transcript
  return(list(gtf_dic,gtf_chm13_to_ensembl_dic))
}
convert_df<-function(opt,gtf_dic){
  df=fread(opt$infile,header=FALSE)
  setnames(df,c("build","sample","enst","count","uniquely_mapped","multi_mapped"))
  print("read")
  df$ensg=gtf_dic[df$enst]
  print("dictionary done")
  print(head(df))
  df_summarised=df%>%
                    group_by(ensg,build,sample)%>%
                    summarise(sum_count=sum(count),
                                sum_uniquely_mapped=sum(uniquely_mapped),
                                sum_multi_mapped=sum(multi_mapped),
                                sum_num_transcripts=n())%>%
                    mutate(prop_uniquely_mapped=sum_uniquely_mapped/sum_count,
                         prop_multi_mapped=sum_multi_mapped/sum_count)
  return(df_summarised)
}
convert_df_chm13<-function(opt,gtf_dics){
  df=fread(opt$infile,header=FALSE)
  setnames(df,c("build","sample","chm13_transcript","count","uniquely_mapped","multi_mapped"))
  chm13_enst_to_ensg=gtf_dics[[1]]
  ensembl_chm13_to_ensg=gtf_dics[[2]]
  df$chm13_gene=chm13_enst_to_ensg[df$chm13_transcript]
  print(head(df))
  df_summarised=df%>%
    group_by(chm13_gene,build,sample)%>%
    summarise(sum_count=sum(count),
              sum_uniquely_mapped=sum(uniquely_mapped),
              sum_multi_mapped=sum(multi_mapped),
              sum_num_transcripts=n())%>%
     mutate(prop_uniquely_mapped=sum_uniquely_mapped/sum_count,
                         prop_multi_mapped=sum_multi_mapped/sum_count)
  df_summarised$ensg=ensembl_chm13_to_ensg[df_summarised$chm13_gene]
  
  return(df_summarised)
}
get_mean<-function(df){
  # df_mean=df%>%group_by(chm13_gene,build)%>%
  #         summarise(mean_uniquely_mapped=mean(sum_uniquely_mapped),
  #                   mean_multi_mapped=mean(sum_multi_mapped),
  #                   mean_num_transcripts=mean(sum_num_transcripts))
  df_mean=df%>%group_by(chm13_gene,build)%>%
          summarise(mean_uniquely_mapped=mean(prop_uniquely_mapped),
                    mean_multi_mapped=mean(prop_multi_mapped),
                    mean_num_counts=mean(sum_count))
  return(df_mean)
}
get_mean_chm13<-function(df){
  # df_mean=df%>%group_by(chm13_gene,build)%>%
  #         summarise(mean_uniquely_mapped=mean(sum_uniquely_mapped),
  #                   mean_multi_mapped=mean(sum_multi_mapped),
  #                   mean_num_transcripts=mean(sum_num_transcripts))
  df_mean=df%>%group_by(chm13_gene,build)%>%
          summarise(mean_uniquely_mapped=mean(prop_uniquely_mapped),
                    mean_multi_mapped=mean(prop_multi_mapped),
                    mean_num_counts=mean(sum_count))
  return(df_mean)
}
main<-function(){
  opt=parse_my_args()
  print(opt)
  if(opt$build=="chm13"){
    print("chm13")
    gtf_dics=get_gtf_chm13(opt)
    print("gtf complete")
    converted_df=convert_df_chm13(opt,gtf_dics)
    print("converted_df complete")
    converted_df_mean=get_mean_chm13(converted_df)
    print("converted_df_mean complete")

  }else{
    print("not chm13")
    gtf_dic=get_gtf(opt)
    print("gtf complete")
    converted_df=convert_df(opt,gtf_dic)
    print("converted_df complete")
    converted_df_mean=get_mean(converted_df)
    print("converted_df_mean complete")

  }
  print("writing converted_df")
  fwrite(converted_df,opt$outfile,quote=FALSE)
  print("writing converted_df_mean")
  fwrite(converted_df_mean,opt$outfile_mean,quote=FALSE)
}

main()




