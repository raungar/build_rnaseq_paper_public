
library(data.table)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

# collate TPM values for all samples of given tissue

# variables
Tissue="Blood"
tissue =  str_to_lower(Tissue)
builds = c("hg19", "hg38") #chm13

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/01_prep_metadata/"
datadir="/oak/stanford/groups/smontgom/shared/UDN/"
metadataFile = paste0(inputdir, "metadata2017to2022_processed_samples.txt")

savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/01_TPM_collated/"
if (!dir.exists(savedir)){
  dir.create(savedir)
} else {
  print("Save directory already exists")
}


# get metadata ----
metadata = read_tsv(metadataFile)
metadata_filt = metadata %>% filter(TISSUE==Tissue)
samples = metadata %>% filter(TISSUE==Tissue) %>% pull(SAMPLE) %>% unique()
n_samples = length(samples)

# get TPM for all samples ----

for (sample in samples){
  
  print(sample)
  
  for (build in builds){
    BUILD = str_to_upper(build)
    # get data
    type="genes"
    rsem_dir = paste0(datadir, "/Preprocessing",BUILD,"Primary/RSEM/")
    suffix=paste0(".sorted_dedup_filtered_rsem_",build,".",type,".results")
    file = paste0(rsem_dir, sample, suffix)
    
    # read and reformat data
    df = read_table(file) %>% 
      # create clean gene name that drops version but retains PAR annotations
      mutate(gene_id = str_remove(gene_id, "_\\d+")) %>%
      mutate(gene = str_remove(gene_id, "\\.\\d+")) %>% 
      mutate(build = paste0(!!build)) %>% 
      mutate(sample = !!sample) %>% 
      mutate(logTPM = log2(TPM + 0.05)) %>% 
      dplyr::select(gene, build, sample, TPM, logTPM)
    # colnames(df)[2] = paste0("tpm_hg",build)
    
    # merge to count matrix
    if (exists("combined")) {
      combined = rbind(combined, df)
    } else {
      combined = df
    }
    
    rm(df)
  }
  
}

# get number of non-expressed genes
n_genes_dropped = combined %>% 
  group_by(gene) %>% 
  filter(sum(TPM > 0.1) <= 5) %>% 
  pull(gene) %>% 
  unique() %>% 
  length()
n_genes_dropped # 28,891 genes with at least 0.1 TPM in 5 or fewer samples

# filter to genes expressed in at least 5 individuals in either build
expressed = combined %>% group_by(gene) %>% filter( sum(TPM > 0.1) > 0.3 * n_samples )
n_genes_exp = expressed %>% pull(gene) %>% unique() %>% length()
n_genes_exp # 35,418 genes with at least 0.1TPM in more than 5 samples

write_tsv(combined, paste0(savedir, tissue, "_tpm_collated.txt"))
write_tsv(expressed, paste0(savedir, tissue, "_tpm_collated_expressed.txt"))



