library(EnsDb.Hsapiens.v86)
library(magrittr)
library(tidyverse)
library(argparser)

# Create a parser
p <- arg_parser("Rank phenotype prioritized genes from solved cases based on expression and splicing zscores")

# Add command line arguments
p <- add_argument(p, "--solved_case_tsv", help="tsv file of solved cases with solved gene", type="character")
p <- add_argument(p, "--ezscore_hg19_blood", help="sample,gene,zscore table for hg19 blood", type="character")
p <- add_argument(p, "--ezscore_hg19_fibro", help="sample,gene,zscore table for hg19 fibroblast", type="character")
p <- add_argument(p, "--ezscore_hg38_blood",help="sample,gene,zscore table for hg38 blood", type="character")
p <- add_argument(p, "--ezscore_hg38_fibro",help="sample,gene,zscore table for hg38 fibroblast", type="character")
p <- add_argument(p,"--ezscore_t2t_blood",help="sample,gene,zscore table for chm13 T2T blood", type="character")
p <- add_argument(p, "--ezscore_t2t_fibro",help="sample,gene,zscore table for chm13 T2T fibroblast", type="character")
p <- add_argument(p, "--szscore_hg19",help="splicing,gene,zscore,tissue table for hg19", type="character")
p <- add_argument(p, "--szscore_hg38",help="splicing,gene,zscore,tissue table for hg38", type="character")
p <- add_argument(p, "--szscore_t2t",help="splicing,gene,zscore,tissue table for t2t", type="character")
p <- add_argument(p, "--gene_list", help="individual,gene table for top 250 phenotype prioritized genes per patient", type="character")

argv <- parse_args(p)


solved_cases <- read_tsv(argv$solved_case_tsv)
solved_cases <- solved_cases[!is.na(solved_cases$GENE),]
zscore <- rbind(
    read_tsv(argv$ezscore_hg19_blood) %>% mutate(BUILD="HG19", Tissue="Blood"),
    read_tsv(argv$ezscore_hg19_fibro) %>% mutate(BUILD="HG19", Tissue="Fibroblast"),
    read_tsv(argv$ezscore_hg38_blood) %>% mutate(BUILD="HG38", Tissue="Blood"),
    read_tsv(argv$ezscore_hg38_fibro) %>% mutate(BUILD="HG38", Tissue="Fibroblast"),
    read_tsv(argv$ezscore_t2t_blood) %>% mutate(BUILD="CHM13", Tissue="Blood"),
    read_tsv(argv$ezscore_t2t_fibro) %>% mutate(BUILD="CHM13", Tissue="Fibroblast")
)
zscore$gene <- gsub(zscore$gene, pattern=c("(ENSG\\d+)\\..*"), replacement="\\1")

splice_zscore <- rbind(
    read_tsv(argv$szscore_hg19, col_names=c("chrom", "start", "end", "cluster_id", "sample_id", "splice_z", "trancsript", "gene", "val")) %>% select(sample_id, cluster_id,gene,splice_z) %>% mutate(BUILD="HG19"),
    read_tsv(argv$szscore_hg38, col_names=c("chrom", "start", "end", "cluster_id", "sample_id", "splice_z", "junc", "col8","col9","splicesite_type","col11", "col12","col13","col14","col15","col16","col17","gene", "col19")) %>% select(sample_id,cluster_id,gene,splice_z) %>% mutate(BUILD="HG38"),
    read_tsv(argv$szscore_t2t, col_names=c("chrom","start","end","cluster_id","sample_id","splice_z","transcript","gene","CHM13_gene", "gene_again","col11")) %>% select(sample_id,cluster_id,gene,splice_z) %>% mutate(BUILD="CHM13")
)
splice_zscore$gene <- gsub(splice_zscore$gene, pattern=c("(ENSG\\d+)\\..*"), replacement="\\1")

gene_lists <- read_tsv(argv$gene_list, col_names=c("Rank", "GENE", "ID", "Score", "Status", "UDNID"),skip=1) %>% dplyr::filter(UDNID %in% solved_cases$UDNID)
gene_symbol <- ensembldb::select(EnsDb.Hsapiens.v86, keys=unique(gene_lists$GENE), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
gene_symbol <- gene_symbol[startsWith(gene_symbol$GENEID, "ENSG"),]
gene_lists %<>% left_join(gene_symbol,  by=c("GENE"="SYMBOL"), multiple="all")
solved_cases <- solved_cases %>% dplyr::filter(UDNID %in% gene_lists$UDNID)
solved_cases$SOLVED_GENE = TRUE
RD_to_UDNID <- solved_cases %>% dplyr::select(RDID, UDNID)
RD_to_UDNID
gene_lists %<>% left_join(RD_to_UDNID, multiple="all")
solved_cases %<>% left_join(gene_lists)
solved_gene_not_prioritized <- solved_cases %>% dplyr::filter(is.na(Rank))

solved_prioritized_gene_df <- distinct(data.frame(UDNID=c(gene_lists$UDNID, solved_gene_not_prioritized$UDNID),
                                 RDID=c(gene_lists$RDID, solved_gene_not_prioritized$RDID),
                                 GENE=c(gene_lists$GENE, solved_gene_not_prioritized$GENE),
                                 GENEID=c(gene_lists$GENEID, solved_gene_not_prioritized$ENSG)))
solved_prioritized_gene_df %<>% left_join(dplyr::select(solved_cases, RDID, GENE, SOLVED_GENE), multiple="all")
solved_prioritized_gene_df$SOLVED_GENE[is.na(solved_prioritized_gene_df$SOLVED_GENE)] <- FALSE
solved_prioritized_gene_df <- solved_prioritized_gene_df[!is.na(solved_prioritized_gene_df$RDID),]
solved_prioritized_gene_df %<>% distinct()
colnames(solved_prioritized_gene_df) <- c("UDNID", "sample_id", "GENE_SYMBOL", "gene", "SOLVED_GENE")

write.table(solved_prioritized_gene_df, file="Solved_cases.RNAseq_samples.pheno_derived_genelist.txt", sep="\t", row.names=F, col.names =T, quote=F)
write.table(prioritized_gene_df, file="Affected_proband.RNAseq_samples.pheno_derived_genelist.txt", sep="\t", row.names=F, col.names =T, quote=F)

solved_prioritized_gene_df_zscores <- solved_prioritized_gene_df %>% left_join(zscore, multiple="all") %>% left_join(splice_zscore)
solved_prioritized_gene_df_zscores <- solved_prioritized_gene_df_zscores[!is.na(solved_prioritized_gene_df_zscores$zscore),]
solved_prioritized_gene_df_zscores %<>% group_by(sample_id, BUILD) %>% mutate(rank_all=rank(rank_all), rank_underexp=rank(rank_underexp), rank_overexp=rank(rank_overexp), rank_splice=rank(-splice_z,ties.method="first"))
solved_prioritized_gene_df_zscores %>% dplyr::filter(SOLVED_GENE) %>% print(n=150)
solved_prioritized_gene_df_zscores %>% dplyr::filter(SOLVED_GENE) %>% print(n=150)
solved_prioritized_gene_df_zscores %>% group_by(sample_id, GENE_SYMBOL) %>% summarize(var=var(rank_splice),n=n(), not_na=!any(is.na(splice_z))) %>% filter(not_na) %>% arrange(desc(var))


table(solved_prioritized_gene_df_zscores$sample_id, solved_prioritized_gene_df_zscores$BUILD)
write.table(solved_prioritized_gene_df_zscores.cleaned, file="Solved_cases.phenotype_derived_genelist.transcriptome_prioritized.hg19_vs_hg38_vs_CHM13.txt", quote=F, row.names=F, col.names =T, sep="\t")

write.table(solved_prioritized_gene_df_zscores)

solved_prioritized_gene_df_zscores.cleaned %>% filter(SOLVED_GENE) %>% pull(sample_id) %>% table %>% length #44 samples
solved_prioritized_gene_df_zscores.cleaned %>% filter(SOLVED_GENE) %>% pull(GENE_SYMBOL) %>% table %>% length #36 diagnostic genes
solved_prioritized_gene_df_zscores.cleaned %>% filter(SOLVED_GENE) %>% pull(UDNID) %>% table %>% length #36 cases
solved_prioritized_gene_df_zscores.cleaned %>% filter(SOLVED_GENE) %>%
    select(sample_id, gene, zscore, BUILD, Tissue) %>% pivot_wider(values_from=zscore, names_from=BUILD) %>%
    ggplot(aes(HG19, HG38)) + geom_point() + theme_classic() + geom_abline(intercept=0, slope=1, color="blue")
ggsave('solved_gene.HG19_vs_hg38.zscore.pdf') 
solved_prioritized_gene_df_zscores %>% filter(SOLVED_GENE) %>%
    select(sample_id, gene, zscore, BUILD, Tissue) %>% pivot_wider(values_from=zscore, names_from=BUILD) %>%
    ggplot(aes(HG38, CHM13)) + geom_point() + theme_classic() + geom_abline(intercept=0, slope=1, color="blue")
ggsave('solved_gene.hg38_vs_chm13.zscore.pdf')
solved_prioritized_gene_df_zscores %>% filter(SOLVED_GENE) %>%
    select(sample_id, gene, rank_underexp, BUILD, Tissue) %>% pivot_wider(values_from=rank_underexp, names_from=BUILD) %>%
    ggplot(aes(HG19, HG38)) + geom_point() + theme_classic() + geom_abline(intercept=0, slope=1, color="blue")
ggsave('solved_gene.HG19_vs_hg38.rank_underexp.pdf')
solved_prioritized_gene_df_zscores %>% filter(SOLVED_GENE) %>%
    select(sample_id, gene, rank_underexp, BUILD, Tissue) %>% pivot_wider(values_from=rank_underexp, names_from=BUILD) %>%
    ggplot(aes(HG38, CHM13)) + geom_point() + theme_classic() + geom_abline(intercept=0, slope=1, color="blue")
ggsave('solved_gene.HG19_vs_chm13.rank_underexp.pdf')

solved_prioritized_gene_df_zscores %>% filter(SOLVED_GENE) %>%
    select(sample_id, gene, rank_all, BUILD, Tissue) %>% pivot_wider(values_from=rank_all, names_from=BUILD) %>%
    ggplot(aes(HG19, HG38)) + geom_point() + theme_classic() + geom_abline(intercept=0, slope=1, color="blue")
ggsave('solved_gene.HG19_vs_hg38.rank_all.pdf')
solved_prioritized_gene_df_zscores %>% filter(SOLVED_GENE) %>%
    select(sample_id, gene, rank_all, BUILD, Tissue) %>% pivot_wider(values_from=rank_all, names_from=BUILD) %>%
    ggplot(aes(HG38, CHM13)) + geom_point() + theme_classic() + geom_abline(intercept=0, slope=1, color="blue")
ggsave('solved_gene.HG19_vs_chm13.rank_all.pdf')
