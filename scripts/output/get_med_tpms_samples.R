library(data.table)
library(tidyverse)


builds='hg38'
tissues=c("Blood","Bone_marrow","Cardiac_tissue","Cerebellum","Embryonic_bodies","Fibroblast","iPSC_NPC","iPS","Muscle","PBMC","Skin")
dir_prefix='/Volumes/groups/smontgom/shared/UDN/Output/hg38/eOutliers'
hgnc_file<-"/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/hgnc_ensg_convert_complete_annotated.txt"


all_tpms=data.frame()
for(build in builds){
  suffix=paste0('.',build,'.dedupOptical_minMQ255_rsem.genes.results')
  for(tissue in tissues){
    this_file=paste0(dir_prefix,'/',tissue,suffix)
    tpms=fread(this_file)
    colnames(tpms)<-c('sample','gene','transcript','tpm')
    tpms_summarized=tpms%>%group_by(gene)%>%summarise(median_tpm=median(tpm),mean_tpm=mean(tpm))%>%mutate(tissue=tissue,num_samples=length(unique(tpms$sample)))
    all_tpms=rbind(all_tpms,tpms_summarized)
  }
}

all_tpms$gene=gsub('\\..*','',all_tpms$gene)
hgnc_convert=fread(hgnc_file)

all_tpms_anno=merge(all_tpms,hgnc_convert[,c('ensg','hgnc','gene_biotype','chrom')],by.x='gene',by.y='ensg',all.x=T)
fwrite(all_tpms_anno%>%arrange(gene),'/Volumes/groups/smontgom/shared/UDN/OutputForGCsHG38Primary/median_mean_TPMs_all_tissues.txt')
