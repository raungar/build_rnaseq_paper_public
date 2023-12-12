
library(data.table)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

# correlate gene expression across individuals -  
# which genes are least concordant across dataset between builds?

# get metadata ----
Tissue="Blood"
tissue =  str_to_lower(Tissue)
builds = c("hg19", "hg38") #chm13

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/LIMMA_DREAM/output/01_prep_metadata/"
datadir="/oak/stanford/groups/smontgom/shared/UDN/"
metadataFile = paste0(inputdir, "metadata2017to2022_processed_samples.txt")

metadata = read_tsv(metadataFile)
metadata_filt = metadata %>% filter(TISSUE==Tissue)
samples = metadata %>% filter(TISSUE==Tissue) %>% pull(SAMPLE) %>% unique()

# get TPM for all samples ----
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01_TPM_all_samples/"

## BLOOD ====
for (sample in samples_blood){
  
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
expressed = combined %>% group_by(gene) %>% filter(sum(TPM > 0.1) > 5)
n_genes_exp = expressed %>% pull(gene) %>% unique() %>% length()
n_genes_exp # 35,418 genes with at least 0.1TPM in more than 5 samples

write_tsv(combined, paste0(savedir, tissue, "_tpm_collated_jul2022.txt"))
write_tsv(expressed, paste0(savedir, tissue, "_tpm_collated_jul2022_expressed.txt"))

rm(combined)


## FIBROBLAST ====
datadir="/oak/stanford/groups/smontgom/shared/UDN/"
for (sample in samples_fibro){
  
  print(sample)
  
  for (build in c("19","38")){
    # get data
    type="genes"
    rsem_dir = paste0(datadir, "/PreprocessingHG",build,"Primary/RSEM/")
    suffix=paste0(".sorted_dedup_filtered_rsem_hg",build,".",type,".results")
    file = paste0(rsem_dir, sample, suffix)
    
    # read and reformat data
    df = read_table(file) %>% 
      # create clean gene name that drops version but retains PAR annotations
      mutate(gene_id = str_remove(gene_id, "_\\d+")) %>%
      mutate(gene = str_remove(gene_id, "\\.\\d+")) %>% 
      mutate(build = paste0("hg", !!build)) %>% 
      mutate(sample = !!sample) %>% 
      mutate(logTPM = log2(TPM + 0.05)) %>% 
      dplyr::select(gene, build, sample, TPM, logTPM)
    # colnames(df)[2] = paste0("tpm_hg",build)
    
    # merge to count matrix
    if (exists("fibro")) {
      fibro = rbind(fibro, df)
    } else {
      fibro = df
    }
    
    rm(df)
  }
  
}

# get number of non-expressed genes
n_genes_dropped = fibro %>% 
  group_by(gene) %>% 
  filter(sum(TPM > 0.1) <= 5) %>% 
  pull(gene) %>% 
  unique() %>% 
  length()
n_genes_dropped # 35,887 genes with at least 0.1 TPM in 5 or fewer samples

# filter to genes expressed in at least 5 individuals in either build
fibro_exp = fibro %>% group_by(gene) %>% filter(sum(TPM > 0.1) > 5)
n_genes_exp = fibro_exp %>% pull(gene) %>% unique() %>% length()
n_genes_exp # 28,422 genes with at least 0.1TPM in more than 5 samples

write_tsv(fibro, paste0(savedir,"output/fibro_tpm_collated_jul2022.txt"))
write_tsv(fibro_exp, paste0(savedir,"output/fibro_tpm_collated_jul2022_expressed.txt"))

rm(fibro)


# Correlations ----
inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"

file = paste0(inputdir, "blood_tpm_collated_jul2022_expressed.txt")
blood_exp = fread(file, header=T, sep='\t')

file = paste0(inputdir, "fibro_tpm_collated_jul2022_expressed.txt")
fibro_exp = fread(file, header=T, sep='\t')

rm(file)

## raw TPM ====

### BLOOD ----
blood_tpm = blood_exp %>% 
  dplyr::select(gene, sample, build, TPM) %>% 
  spread(key = build, value = TPM) %>% 
  rename(hg19=`19`, hg38=`38`)

blood_tpm_cor = blood_tpm %>% 
  group_by(gene) %>% 
  summarize(pearson = cor(hg19, hg38, method = 'pearson'),
            spearman = cor(hg19, hg38, method = 'spearman'),
            covar = var(hg19, hg38, na.rm=T), 
            n_TPMgt0_hg38 = sum(hg38>0), 
            n_TPMgt0_hg19 = sum(hg19>0),
            TPM_median_hg38 = median(hg38),
            TPM_median_hg19 = median(hg19),
            logFC_median = median(log(hg38+0.05)-log(hg19+0.05)),
            logFC_mean = mean(log(hg38+0.05)-log(hg19+0.05)),
            logFC_sd = stats::sd(log(hg38+0.05)-log(hg19+0.05))) %>% 
  arrange(abs(pearson))

## save results
write_tsv(blood_tpm, paste0(savedir, "blood_tpm_collated_jul2022_expressed.TPM.txt"))
write_tsv(blood_tpm_cor, paste0(savedir, "blood_tpm_collated_jul2022_expressed.TPM_correlations.txt"))

### FIBROBLAST ----
fibro_tpm = fibro_exp %>% 
  dplyr::select(gene, sample, build, TPM) %>% 
  spread(key = build, value = TPM) 

fibro_tpm_cor = fibro_tpm %>% 
  group_by(gene) %>% 
  summarize(pearson = cor(hg19, hg38, method = 'pearson'),
            spearman = cor(hg19, hg38, method = 'spearman'),
            covar = var(hg19, hg38, na.rm=T), 
            n_TPMgt0_hg38 = sum(hg38>0), 
            n_TPMgt0_hg19 = sum(hg19>0),
            TPM_median_hg38 = median(hg38),
            TPM_median_hg19 = median(hg19),
            logFC_median = median(log(hg38+0.05)-log(hg19+0.05)),
            logFC_mean = mean(log(hg38+0.05)-log(hg19+0.05)),
            logFC_sd = stats::sd(log(hg38+0.05)-log(hg19+0.05))) %>% 
  arrange(abs(pearson))

## save results
write_tsv(fibro_tpm, paste0(savedir, "fibro_tpm_collated_jul2022_expressed.TPM.txt"))
write_tsv(fibro_tpm_cor, paste0(savedir, "fibro_tpm_collated_jul2022_expressed.TPM_correlations.txt"))


## logTPM ====

### BLOOD ----
blood_logtpm = blood_exp %>% 
  dplyr::select(gene, sample, build, logTPM) %>% 
  spread(key = build, value = logTPM)%>% 
  rename(hg19=`19`, hg38=`38`)
  

# rm(blood_exp, blood_tpm)

blood_logtpm_cor = blood_logtpm %>% 
  group_by(gene) %>% 
  summarize(pearson = cor(hg19, hg38, method = 'pearson'),
            spearman = cor(hg19, hg38, method = 'spearman'),
            covar = var(hg19, hg38, na.rm=T), 
            n_TPMgt0_hg38 = sum(hg38>log2(0 + 0.05) & !is.na(hg38)),
            n_TPMgt0_hg19 = sum(hg19>log2(0 + 0.05) & !is.na(hg19)),
            logTPM_median_hg38 = median(hg38),
            logTPM_median_hg19 = median(hg19),
            logFC_median = median(hg38-hg19),
            logFC_mean = mean(hg38-hg19),
            logFC_sd = stats::sd(hg38-hg19)) %>% 
  arrange(abs(pearson)) %>% 
  # mutate(n_diff = n_TPMgt0_hg38 - n_TPMgt0_hg19) %>% 
  mutate(logTPM_median_diff = logTPM_median_hg38 - logTPM_median_hg19)

## save results
write_tsv(blood_logtpm, paste0(savedir, "blood_tpm_collated_jul2022_expressed.logTPM.txt"))
write_tsv(blood_logtpm_cor, paste0(savedir, "blood_tpm_collated_jul2022_expressed.logTPM_correlations.txt"))


### FIBROBLAST ----
fibro_logtpm = fibro_exp %>% 
  dplyr::select(gene, sample, build, logTPM) %>% 
  spread(key = build, value = logTPM)


# rm(blood_exp, blood_tpm)

fibro_logtpm_cor = fibro_logtpm %>% 
  group_by(gene) %>% 
  summarize(pearson = cor(hg19, hg38, method = 'pearson'),
            spearman = cor(hg19, hg38, method = 'spearman'),
            covar = var(hg19, hg38, na.rm=T), 
            n_TPMgt0_hg38 = sum(hg38>log2(0 + 0.05) & !is.na(hg38)),
            n_TPMgt0_hg19 = sum(hg19>log2(0 + 0.05) & !is.na(hg19)),
            logTPM_median_hg38 = median(hg38),
            logTPM_median_hg19 = median(hg19),
            logFC_median = median(hg38-hg19),
            logFC_mean = mean(hg38-hg19),
            logFC_sd = stats::sd(hg38-hg19)) %>% 
  arrange(abs(pearson)) %>% 
  # mutate(n_diff = n_TPMgt0_hg38 - n_TPMgt0_hg19) %>% 
  mutate(logTPM_median_diff = logTPM_median_hg38 - logTPM_median_hg19)

## save results
write_tsv(fibro_logtpm, paste0(savedir, "fibro_tpm_collated_jul2022_expressed.logTPM.txt"))
write_tsv(fibro_logtpm_cor, paste0(savedir, "fibro_tpm_collated_jul2022_expressed.logTPM_correlations.txt"))


### Sanity: Find samples that are uniquely discordant
blood_logtpm_cor_samples = blood_logtpm %>% 
  group_by(sample) %>% 
  drop_na() %>% 
  summarize(pearson = cor(hg19, hg38, method = 'pearson'),
            spearman = cor(hg19, hg38, method = 'spearman'),
            covar = var(hg19, hg38, na.rm=T),
            n_genes = n()) %>% 
  arrange(abs(spearman)) 

summary(blood_logtpm_cor_samples$pearson) # min: 0.9907  max: 0.9953
summary(blood_logtpm_cor_samples$spearman) # min: 0.9912  max: 0.9951
rm(blood_logtpm_cor_samples)

fibro_logtpm_cor_samples = fibro_logtpm %>% 
  group_by(sample) %>% 
  drop_na() %>% 
  summarize(pearson = cor(hg19, hg38, method = 'pearson'),
            spearman = cor(hg19, hg38, method = 'spearman'),
            covar = var(hg19, hg38, na.rm=T),
            n_genes = n()) %>% 
  arrange(abs(spearman)) 

summary(fibro_logtpm_cor_samples$pearson) # min: 0.9930  max: 0.9956
summary(fibro_logtpm_cor_samples$spearman) # min: 0.9933  max: 0.9958
rm(fibro_logtpm_cor_samples)


# visualize results ----

## BLOOD ----

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/output/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"

tissue="blood"

logtpm_file = paste0(inputdir, tissue, "_tpm_collated_jul2022_expressed.logTPM.txt")
logtpm_cor_file = paste0(inputdir, tissue, "_tpm_collated_jul2022_expressed.logTPM_correlations.txt")
logtpm = read_tsv(logtpm_file)
logtpm_cor = read_tsv(logtpm_cor_file)

n_genes = logtpm_cor %>% filter(!is.na(pearson)) %>% pull(gene) %>% unique() %>% length()
n_samples = logtpm %>% pull(sample) %>% unique() %>% length()

### check correlations ====
logtpm_cor %>% 
  filter(pearson < 0.5 & spearman < 0.5) %>% 
  nrow() # 157 genes with low concordance in TPM between hg19 and hg38


### TPM scatter ----

# logtpm_scatter = fibro_logtpm %>%
#   ggplot() +
#   geom_point(aes(x = hg19, y = hg38), alpha = 0.3, size = 1, shape=20) +
#   theme_bw(base_size = 14)
# # logtpm_scatter
# # 
# logtpm_scatter + facet_wrap (~sample)

### Distributions ----


plot_corr_distr = logtpm_cor %>% 
  gather(key = "method", value = "correlation", pearson, spearman) %>% 
  ggplot(aes(x = correlation, fill = method)) +
  geom_histogram() +
  scale_fill_viridis_d(option = "E", guide="none") +
  facet_wrap(~ method) +
  theme_classic(base_size = 20) +
  ggtitle("Correlation coefficient distributions from hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " ",tissue," RNA-seq samples from the UDN"))

plot_corr_distr
plot_corr_distr + facet_wrap(~ method, scales="free") + scale_y_continuous(trans="log10", labels=scales::comma)

ggsave(plot=plot_corr_distr, filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.histogram_linear.png"))
ggsave(plot=plot_corr_distr + scale_y_continuous(trans="log10", labels=scales::comma), 
       filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.histogram_log.png"))


### Coef scatter ----

plot_corr_compare = logtpm_cor %>%
  ggplot(aes(x = spearman, y = pearson, size = abs(logFC_median), color = logFC_median)) +
  geom_abline(slope = 1, color = "grey17", linetype = "dashed", alpha = 0.2) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(name  = 'median logFC\nlog(TPM_hg38) - log(TPM_hg19)') +
  scale_size_continuous(name = "absolute median logFC") +
  theme_classic(base_size = 18) +
  ggtitle("Pearson vs. Spearman correlation between hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " ",tissue," RNA-seq samples from the UDN"))

plot_corr_compare

ggsave(plot=plot_corr_compare, filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.spearman_vs_pearson.png"))


# get discordant genes
logtpm_cor %>% filter(spearman < 0.5 & pearson < 0.5) %>% nrow() #157

### Pearson vs logFCmedian ----

plot_corr_logFC = logtpm_cor %>% 
  ggplot(aes(x = logFC_median, y = -log10(pearson), size = abs(logFC_median), color = logFC_median)) +
  geom_hline(yintercept = -log10(0.5), linetype="dashed", color="red") +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(name  = 'median logFC\nlog(TPM_hg38) - log(TPM_hg19)') +
  scale_size_continuous(name = "absolute median logFC") +
  theme_classic(base_size = 16) +
  ggtitle("Pearson vs. Median LogFC between hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " ",tissue," RNA-seq samples from the UDN"))

ggsave(plot=plot_corr_logFC, filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.pearson_vs_medLogFC.png"))


## FIBROBLASTS ----

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"

logtpm_file = paste0(inputdir, "tpm_collated_jul2022_expressed.logTPM.txt")
logtpm_cor_file = paste0(inputdir, "tpm_collated_jul2022_expressed.logTPM_correlations.txt")
logtpm = read_tsv(logtpm_file)
logtpm_cor = read_tsv(logtpm_cor_file)


### check correlations ====
blood_logtpm_cor %>% 
  filter(pearson < 0.5 & spearman < 0.5) %>% 
  nrow() # 157 genes with low concordance in TPM between hg19 and hg38


logtpm_cor %>% 
  filter(pearson < 0.5 & spearman < 0.5) %>% 
  nrow() # 79 genes with low concordance in TPM between hg19 and hg38



### TPM scatter ----

logtpm_scatter = logtpm %>%
  ggplot() +
  geom_point(aes(x = hg19, y = hg38), alpha = 0.3, size = 1, shape=20) +
  theme_bw(base_size = 14)
# logtpm_scatter
# 
logtpm_scatter + facet_wrap (~sample)

### Distributions ----

n_genes = blood_logtpm_cor %>% filter(!is.na(pearson)) %>% pull(gene) %>% unique() %>% length()
n_samples = blood_logtpm %>% pull(sample) %>% unique() %>% length()

plot_corr_distr = blood_logtpm_cor %>% 
  gather(key = "method", value = "correlation", pearson, spearman) %>% 
  ggplot(aes(x = correlation, fill = method)) +
  geom_histogram() +
  scale_fill_viridis_d(option = "E", guide="none") +
  facet_wrap(~ method) +
  theme_classic(base_size = 20) +
  ggtitle("Correlation coefficient distributions from hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " ",tissue," RNA-seq samples from the UDN"))

plot_corr_distr
plot_corr_distr + facet_wrap(~ method, scales="free") + scale_y_continuous(trans="log10", labels=scales::comma)

ggsave(plot=plot_corr_distr, filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.histogram_linear.png"))
ggsave(plot=plot_corr_distr + scale_y_continuous(trans="log10", labels=scales::comma), 
       filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.histogram_log.png"))


### Coef scatter ----

plot_corr_compare = blood_logtpm_cor %>%
  ggplot(aes(x = spearman, y = pearson, size = abs(logFC_median), color = logFC_median)) +
  geom_abline(slope = 1, color = "grey17", linetype = "dashed", alpha = 0.2) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(name  = 'median logFC\nlog(TPM_hg38) - log(TPM_hg19)') +
  scale_size_continuous(name = "absolute median logFC") +
  theme_classic(base_size = 18) +
  ggtitle("Pearson vs. Spearman correlation between hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " ",tissue," RNA-seq samples from the UDN"))

plot_corr_compare

ggsave(plot=plot_corr_compare, filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.spearman_vs_pearson.png"))


# get discordant genes
blood_logtpm_cor %>% filter(spearman < 0.5 & pearson < 0.5) %>% nrow() #157

### Pearson vs logFCmedian ----

plot_corr_logFC = blood_logtpm_cor %>% 
  ggplot(aes(x = logFC_median, y = -log10(pearson), size = abs(logFC_median), color = logFC_median)) +
  geom_hline(yintercept = -log10(0.5), linetype="dashed", color="red") +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(name  = 'median logFC\nlog(TPM_hg38) - log(TPM_hg19)') +
  scale_size_continuous(name = "absolute median logFC") +
  theme_classic(base_size = 16) +
  ggtitle("Pearson vs. Median LogFC between hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " ",tissue," RNA-seq samples from the UDN"))

ggsave(plot=plot_corr_logFC, filename=paste0(savedir, "figures/",tissue,"_tpm_collated_jul2022.logTPM_correlations.pearson_vs_medLogFC.png"))


# Annotate ----

library("biomaRt")

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
#blood_logtpm_file = paste0(inputdir, "blood_tpm_collated_jul2022_expressed.logTPM.txt")
blood_logtpm_cor_file = paste0(inputdir, "blood_tpm_collated_jul2022_expressed.logTPM_correlations.txt")
#blood_logtpm = read_tsv(blood_logtpm_file)
blood_logtpm_cor = read_tsv(blood_logtpm_cor_file)


## Gene Names ====

mart <- useEnsembl(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl")
mybm=getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
           mart = mart)

# create hgnc- > ensembl dictionary
hgnc_to_ensg <- mybm$ensembl_gene_id
names(hgnc_to_ensg) <- mybm$hgnc_symbol
hgnc_to_ensg["CBSL"]<-hgnc_to_ensg["CBS"]
hgnc_to_ensg["SMIM11A"]<-hgnc_to_ensg["SMIM11"]
hgnc_to_ensg["SMIM11B"]<-hgnc_to_ensg["SMIM11"]

# create ensembl > hgnc dictionary
ensg_to_hgnc <- names(hgnc_to_ensg)
names(ensg_to_hgnc) = hgnc_to_ensg


# annotate
blood_logtpm_cor_annot = blood_logtpm_cor %>% 
  mutate(gene_name = ensg_to_hgnc[gene])

## Genes with different versions ====
# read and reformat discrepant gene version data
discrep_vers = read_tsv(paste0(inputdir, "output/02_prep_countdata/Blood_different_gene_versions_hg38_hg19.txt"))
colnames(discrep_vers) = c('gene', 'gene_version', 'hg19', 'hg38')
discrep_vers = discrep_vers %>% 
  gather(key = 'build', value = 'count', hg19, hg38) %>% 
  drop_na() %>% 
  dplyr::select(-count) %>% 
  spread(key = build, value = gene_version)
head(discrep_vers)

blood_logtpm_cor_annot = blood_logtpm_cor_annot %>% 
  mutate(diff_version = ifelse(gene %in% discrep_vers$gene, "Yes", NA)) %>% 
  left_join(discrep_vers %>% 
              rename(hg19_version = hg19, hg38_version = hg38))


## DEGenes ====
inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
degenes = read_tsv(paste0(inputdir, "output/03_diff_expression/Blood_diff_expression.limma_voom_dream.all.txt"))
colnames(degenes) <- c('gene', paste("LIMMA", colnames(degenes)[-1], sep = "_"))

# merge
blood_logtpm_cor_annot = left_join(blood_logtpm_cor_annot, degenes)

# get top discrepant genes
head = blood_logtpm_cor_annot %>% filter( (spearman < 0.5 & pearson < 0.5) | 
                                          (abs(LIMMA_logFC) > 1.5 & LIMMA_adj.P.Val < 0.05) )



## Discrepant Genes ====
ref_dir="/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/"
discreps_genes = fread(paste0(ref_dir, "exome_discrep_enriched_genes_li2021.txt")) %>% 
  mutate(gene = hgnc_to_ensg[gene])

blood_logtpm_cor_annot = blood_logtpm_cor_annot %>% 
  left_join(discreps_genes %>%
              mutate(discrep_var_enriched = "Yes") %>% 
              dplyr::select(gene, discrep_var_enriched), 
            by = c("gene"="gene"))

## OMIM Genes ====
ref_dir="/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/"
pheno2gene="morbidmap.txt"
pheno = fread(paste0(ref_dir, pheno2gene))
colnames(pheno) = c("pheno","gene","mim","cyto")

pheno_clean = pheno %>% 
  extract(pheno, into = c("pheno", "confidence"), "(.*) \\(([0-9]+)\\)$") %>%
  mutate(confidence = as.numeric(confidence)) %>% 
  extract(pheno, into = c("pheno", "pheno_mim"), "(.*), ([0-9]{6})$") %>% 
  mutate(status = ifelse(confidence == 1, "association", 
                         ifelse(confidence == 2, "linakge_mapping", 
                                ifelse(confidence == 3, "known_molec_cause", 
                                       ifelse(confidence == 4, "del_dup_syndrome", NA))))) %>% 
  separate_rows(gene, convert = T, sep=", ") %>% 
  dplyr::select(-mim) %>% 
  mutate(ensg = hgnc_to_ensg[gene]) %>% 
  distinct() %>% 
  # add something to fill down ensg within groups
  group_by(pheno, pheno_mim, confidence, status, cyto) %>% 
  mutate(approved_symbol = ifelse(!is.na(ensg), gene, NA)) %>% 
  summarize(ensg = toString(na.omit(ensg)),
            approved_symbol = toString(na.omit(approved_symbol)),
            hgnc = toString(na.omit(gene))) %>% 
  ungroup() %>% 
  separate_rows(ensg, approved_symbol, convert=T, sep=", ") %>% 
  group_by(ensg, hgnc, approved_symbol) %>% 
  summarize(pheno = toString(unique(na.omit(pheno))),
            pheno_mim  = toString(unique(na.omit(pheno_mim ))),
            confidence = toString(unique(na.omit(confidence))),
            status = toString(unique(na.omit(status)))) %>% 
  ungroup()

blood_logtpm_cor_annot = blood_logtpm_cor_annot %>% 
  left_join(pheno_clean %>% dplyr::select(gene = ensg, 
                                          OMIM_trait = pheno, 
                                          OMIM_pheno_id = pheno_mim, 
                                          OMIM_confidence = confidence, 
                                          OMIM_status = status), 
            by = c("gene"="gene")) %>% 
  mutate(OMIM_gene = ifelse(OMIM_confidence >= 3, "Yes", NA)) 

## ANI ====
ani_dir="~/UDN/UDNShared/Analysis/ReferenceComparison/RefGenomeANI/output/hg19_vs_hg38/"
ani_file=paste0(ani_dir, "hg19_vs_hg38.exonic.ANI.txt")
ani = fread(ani_file) %>% dplyr::select(gene = geneid, ANI)

blood_logtpm_cor_annot = blood_logtpm_cor_annot %>% left_join(ani)

## HG19 - HG38 Known differences

# COMING SOON #

## Save annotations ====
blood_logtpm_cor_annot = blood_logtpm_cor_annot %>% 
  dplyr::select(gene, gene_name, pearson:logTPM_median_diff, ends_with('version'),
                starts_with('LIMMA'), starts_with('OMIM'), discrep_var_enriched, ANI) %>% 
  distinct()

save_file = paste0(str_remove(blood_logtpm_cor_file, "txt$"), "annotated.txt")
write_tsv(blood_logtpm_cor_annot, save_file)

## Table ====
blood_logtpm_cor_annot %>% 
  filter(pearson < 0.5) %>%
  filter(abs(logFC_median) >= 1.5) %>%
  summarize(total = n(),
            diff_version = sum(!is.na(diff_version)),
            discrep_var_enriched = sum(!is.na(discrep_var_enriched)),
            OMIM_gene = sum(!is.na(OMIM_gene)))

  

# Joint Plots ----
library(ggrepel)

## annotated volcano plot ====
annotated_volcano = blood_logtpm_cor_annot %>% 
  mutate(OMIM_gene = ifelse(is.na(OMIM_gene), "No", "Yes")) %>% 
  {
  ggplot(data = .,
         aes(x = logFC_median, 
             y = -log10(pearson))
         ) +
  geom_hline(yintercept = -log10(0.5), linetype="dashed", color="red") +
      
  geom_point(aes(size = abs(logFC_median), color = ANI,
                 shape = OMIM_gene, alpha = OMIM_gene)) +
      
  geom_label_repel(data = filter(., abs(logFC_median) >= 1.5 & 
                                   OMIM_gene == "Yes" &
                                   abs(pearson) <= 0.5),
                   aes(label = gene_name),
                   point.padding = 0.5, 
                   min.segment.length = 0,
                   nudge_x = -3,
                   nudge_y = 0.2) +
      
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_color_viridis_c(name  = 'Codon Nucleotide\nSimilarity') +
  scale_size_continuous(name = "absolute median logFC") +
      
  theme_classic(base_size = 16) +
      
  ggtitle("Pearson vs. Median LogFC between hg38 and hg19 log(TPM) values", 
          subtitle = paste0("for each of ", n_genes ," genes based on ", n_samples,
                            " blood samples"))
}

annotated_volcano
ggsave(plot=annotated_volcano, filename=paste0(savedir, "blood_tpm_collated_jul2022.logTPM_correlations.pearson_volcano_annot.png"))

## compare to LIMMA FC ====
### sanity check - logFC_median from correlation script and logFC from LIMMA should be well correlated
limma_vs_median = blood_logtpm_cor_annot %>% 
  mutate(`Abs(logFC) > 1.5 in` = ifelse(abs(LIMMA_logFC) >= 1.5 & abs(logFC_median) < 1.5, 
                                        "LIMMA only", 
                                        ifelse(abs(LIMMA_logFC) < 1.5 & abs(logFC_median) >= 1.5,
                                               "Median logFC Only", 
                                               ifelse(abs(LIMMA_logFC) >= 1.5 & abs(logFC_median) > 1.5, "Both", "Neither"))
  )
  ) %>% 
  ggplot(aes(x = logFC_median, y = LIMMA_logFC, color = `Abs(logFC) > 1.5 in`)) +
  geom_abline(slope = 1, linetype="dashed", color="grey67") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("darkolivegreen3", "darkorchid","darkturquoise","grey57")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.8,0.2)) +
  ggtitle("LIMMA-calculated logFC vs. Median logFC")

limma_vs_median
ggsave(plot=limma_vs_median, filename=paste0(savedir, "blood_tpm_collated_jul2022.logTPM_correlations.LIMMAlogFC_vs_medlogFC.png"))

## Compare to LIMMA Sig Genes ====
blood_logtpm_cor_annot %>% 
  filter(LIMMA_adj.P.Val < 0.05) %>%
  filter(abs(LIMMA_logFC) >= 1.5) %>% 
  ggplot(aes(x=pearson)) + geom_histogram()
  
blood_logtpm_cor_annot %>% 
  filter(pearson < 0.5) %>%
  filter(abs(logFC_median) >= 1.5) %>% 
  ggplot(aes(x=LIMMA_logFC)) + geom_histogram()

