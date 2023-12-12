
library(data.table)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

# correlate gene expression across individuals -  
# which genes are least concordant across dataset between builds?

# variables
Tissue="Blood"
tissue =  str_to_lower(Tissue)
builds = c("hg19", "hg38") #chm13

# inputdir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/01_TPM_collated/"
inputdir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/output/06_TPM_pearson/"
if (!dir.exists(savedir)){
  dir.create(savedir)
} else {
  print("Save directory already exists")
}


# Correlations ----

file = paste0(inputdir, tissue, "_tpm_collated_expressed.txt")
expressed = fread(file, header=T, sep='\t')

rm(file)

## raw TPM ====

tpm = expressed %>% 
  dplyr::select(gene, sample, build, TPM) %>% 
  spread(key = build, value = TPM) %>% 
  rename(hg19=`19`, hg38=`38`)

tpm_cor = tpm %>% 
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
write_tsv(tpm, paste0(savedir, tissue, "_tpm_collated_expressed.TPM.txt"))
write_tsv(tpm_cor, paste0(savedir, tissue, "_tpm_collated_expressed.TPM_correlations.txt"))

## logTPM ====

logtpm = expressed %>% 
  dplyr::select(gene, sample, build, logTPM) %>% 
  spread(key = build, value = logTPM) %>% 
  rename(hg19=`19`, hg38=`38`) %>%
  group_by(gene) %>% 
  mutate(allNA = !all(is.na(hg19)|is.na(hg38))) %>%
  filter(allNA) %>%
  dplyr::select(-allNA)


# rm(expressed, blood_tpm)

logtpm_cor = logtpm %>% 
  group_by(gene) %>% 
  summarize(slope = unname(lm(hg19~hg38, na.action=na.omit)$coefficients[2]),
            r2 = summary(lm(hg19~hg38, na.action=na.omit))$r.squared,
            # resid_sd = summary(lm(hg19~hg38, na.action=na.omit))$r.squared,
            pearson = cor(hg19, hg38, method = 'pearson'),
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
write_tsv(logtpm, paste0(savedir, tissue, "_tpm_collated_expressed.logTPM.txt"))
write_tsv(logtpm_cor, paste0(savedir, tissue, "_tpm_collated_expressed.logTPM_correlations.txt"))


# visualize results ----

## BLOOD ----

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/output/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"

tissue="blood"

logtpm_file = paste0(inputdir, tissue, "_tpm_collated_expressed.logTPM.txt")
logtpm_cor_file = paste0(inputdir, tissue, "_tpm_collated_expressed.logTPM_correlations.txt")
logtpm = read_tsv(logtpm_file)
logtpm_cor = read_tsv(logtpm_cor_file)

n_genes = logtpm_cor %>% filter(!is.na(pearson)) %>% pull(gene) %>% unique() %>% length()
n_samples = logtpm %>% pull(sample) %>% unique() %>% length()

### check correlations ====
logtpm_cor %>% 
  filter(r2 < 0.5 ) %>% 
  nrow() # 364 genes with low correlation in TPM between hg19 and hg38

logtpm_cor %>% 
  filter(r2 < 0.5 & (slope < 0.5 | slope > 1.5) ) %>% 
  nrow() # 193 genes with low correlation in TPM between hg19 and hg38

logtpm_cor %>% 
  ggplot(aes(y = slope, x = r2, color = slope-1, alpha = 1-r2)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey56") + 
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey56") + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") + 
  geom_point(size = 2, alpha = 0.5) +
  scale_colour_distiller(
    type = "div",
    palette = 1,
    direction = 1,
    aesthetics = "colour"
  ) +
  scale_alpha_continuous(range = c(0.5, 0.9)) +
  theme_bw(base_size = 14) +
  xlab("Pearson Correlation r^2") + 
  ylab("Slope of line of best fit") +
  ggtitle(paste0("Correlation of hg38 and hg19 logTPM values\n", 
  "for ", n_genes, " across ", n_samples, " ", tissue, " samples"))

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
  dplyr::select(-resid_sd, -pearson) %>% 
  gather(key = "method", value = "correlation", slope, r2, covar) %>% 
  ggplot(aes(x = correlation, fill = method)) +
  geom_histogram() +
  scale_fill_viridis_d(option = "E", guide="none") +
  facet_wrap(~ method, scales = "free") +
  scale_y_continuous(trans="log10") +
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
  ggplot(aes(x = slope, y = pearson, size = abs(logFC_median), color = logFC_median)) +
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



# Annotate ----

library("biomaRt")

inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
savedir="/oak/stanford/groups/smontgom/pgoddard/UDN/UDNShared/Analysis/ReferenceComparisonPrimary/TPMDifferences/PearsonCor/"
#logtpm_file = paste0(inputdir, tissue, "_tpm_collated_jul2022_expressed.logTPM.txt")
logtpm_cor_file = paste0(inputdir, tissue, "_tpm_collated_jul2022_expressed.logTPM_correlations.txt")
#logtpm = read_tsv(blood_logtpm_file)
logtpm_cor = read_tsv(blood_logtpm_cor_file)


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
logtpm_cor_annot = logtpm_cor %>% 
  mutate(gene_name = ensg_to_hgnc[gene])

## Genes with different versions ====
# read and reformat discrepant gene version data
discrep_vers = read_tsv(paste0(inputdir, "../LIMMA_DREAM/02_prep_countdata/Blood_different_gene_versions_hg38_hg19.txt"))
colnames(discrep_vers) = c('gene', 'gene_version', 'hg19', 'hg38')
discrep_vers = discrep_vers %>% 
  gather(key = 'build', value = 'count', hg19, hg38) %>% 
  drop_na() %>% 
  dplyr::select(-count) %>% 
  spread(key = build, value = gene_version)
head(discrep_vers)

logtpm_cor_annot = logtpm_cor_annot %>% 
  mutate(diff_version = ifelse(gene %in% discrep_vers$gene, "Yes", NA)) %>% 
  left_join(discrep_vers %>% 
              rename(hg19_version = hg19, hg38_version = hg38))


## DEGenes ====
inputdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/LIMMA_DREAM/"
degenes = read_tsv(paste0(inputdir, "03_diff_expression/Blood_diff_expression.limma_voom_dream.all.txt"))
colnames(degenes) <- c('gene', paste("LIMMA", colnames(degenes)[-1], sep = "_"))

# merge
logtpm_cor_annot = left_join(logtpm_cor_annot, degenes)

# get top discrepant genes
# head = logtpm_cor_annot %>% filter( (spearman < 0.5 & pearson < 0.5) | 
#                                             (abs(LIMMA_logFC) > 1.5 & LIMMA_adj.P.Val < 0.05) )



## Discrepant Genes ====
ref_dir="/oak/stanford/groups/smontgom/shared/UDN/ReferenceFiles/"
discreps_genes = fread(paste0(ref_dir, "exome_discrep_enriched_genes_li2021.txt")) %>% 
  mutate(gene = hgnc_to_ensg[gene])

logtpm_cor_annot = logtpm_cor_annot %>% 
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

logtpm_cor_annot = logtpm_cor_annot %>% 
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

logtpm_cor_annot = logtpm_cor_annot %>% left_join(ani)

## HG19 - HG38 Known differences

# COMING SOON #

## Save annotations ====
logtpm_cor_annot = logtpm_cor_annot %>% 
  dplyr::select(gene, gene_name, slope:logTPM_median_diff, ends_with('version'),
                starts_with('LIMMA'), starts_with('OMIM'), discrep_var_enriched, ANI) %>% 
  distinct()

save_file = paste0(str_remove(logtpm_cor_file, "txt$"), "annotated.txt")
write_tsv(logtpm_cor_annot, save_file)

## Table ====
logtpm_cor_annot %>% 
  filter(r2 < 0.5) %>%
  filter(abs(logFC_median) >= 1.5) %>%
  summarize(total = n(),
            diff_version = sum(!is.na(diff_version)),
            discrep_var_enriched = sum(!is.na(discrep_var_enriched)),
            OMIM_gene = sum(!is.na(OMIM_gene)))



# Joint Plots ----
library(ggrepel)

## annotated volcano plot ====
annotated_volcano = logtpm_cor_annot %>% 
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
limma_vs_median = logtpm_cor_annot %>% 
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
logtpm_cor_annot %>% 
  filter(LIMMA_adj.P.Val < 0.05) %>%
  filter(abs(LIMMA_logFC) >= 1.5) %>% 
  ggplot(aes(x=pearson)) + geom_histogram()

logtpm_cor_annot %>% 
  filter(pearson < 0.5) %>%
  filter(abs(logFC_median) >= 1.5) %>% 
  ggplot(aes(x=LIMMA_logFC)) + geom_histogram()

