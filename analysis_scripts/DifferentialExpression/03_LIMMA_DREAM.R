# DE analysis ====

# TO DO: update this script to give a useful log

suppressMessages(suppressWarnings(library('variancePartition'))) # DREAM
suppressMessages(suppressWarnings(library('edgeR'))) # LIMMA
suppressMessages(suppressWarnings(library('BiocParallel')))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(optparse)))

`%notin%` <- Negate(`%in%`)

cat(paste(Sys.time(), "Starting analysis\n"))

option_list = list(
  make_option(c( "--build1"), type="character", default=NULL, help="name (ex: hg19)", metavar="character"),
  make_option(c("--build2"), type="character", default=NULL, help="name (ex: hg38)", metavar="character"),
  make_option(c("--build1_infile"), type="character", default=NULL, help="compiled tissue.genes.results file for build 1", metavar="character"),
  make_option(c("--build2_infile"), type="character", default=NULL, help="compiled tissue.genes.results file for build 1", metavar="character"),
  make_option(c("--md"), type="character", default=NULL, help="metadata", metavar="character"),
  make_option(c("--tissue"), type="character", default=NULL, help="tissue", metavar="character"),
  make_option(c("--intermediate_outdir"), type="character", default=NULL, help="where to save counts matrix", metavar="character"),  
  make_option(c("--counts_outfile"), type="character", default=NULL, help="where to save counts matrix", metavar="character"),  
  make_option(c("--fitmm_outfile"), type="character", default=NULL, help="where to save fitmm", metavar="character"),
  make_option(c("--sig_outfile"), type="character", default=NULL, help="where to save sig", metavar="character"),
  make_option(c("--results_outfile"), type="character", default=NULL, help="where to save results", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



build1=opt$build1
build2=opt$build2
builds=c(build1, build2)
build1_infile=opt$build1_infile
build2_infile=opt$build2_infile
output_dir1=opt$intermediate_outdir
results_outfile=opt$results_outfile
my_md=opt$md
tissue=opt$tissue


# tissue="PBMC"
# analysis_dir = "/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/"
# input_dir=paste0(analysis_dir, "output/02_PreppedCountData/")
# output_dir1=paste0(analysis_dir, "output/02_intermediateOutputs/")
# output_dir2=paste0(analysis_dir, "output/03_LIMMA_DiffExp_20221112/")
# dir.create(file.path(output_dir1), showWarnings = FALSE)
# dir.create(file.path(output_dir2), showWarnings = FALSE)
# build1 = "chm13"
# build2 = "hg38"
# builds = c(build1, build2)
# build1ens = ifelse(grepl("chm", build1), paste0(build1, "ensembl"), build2)
# build2ens = ifelse(grepl("chm", build2), paste0(build2, "ensembl"), build2)

# build1_infile = paste0(input_dir, tissue, ".", build1ens , ".dedupOptical_minMQ255_rsem.genes.counts_matrix.txt")
# build2_infile = paste0(input_dir, tissue, ".", build2ens , ".dedupOptical_minMQ255_rsem.genes.counts_matrix.txt")

# my_md="/oak/stanford/groups/smontgom/shared/UDN/Data/Metadata/2017to2022_metadata.tsv"

get_metadata<-function(md_file, builds=c(build1, build2), tissue){
  # check that builds is a pair of builds:
  if (length(builds) != 2){
    cat(paste0("error: requires 2 builds. You have entered ", length(builds)))
    stop()
  }
  metadata_full <- read_tsv(md_file)
  ## reformat ----
  metadata_selected = metadata_full %>% 
    dplyr::select(INDIVIDUAL = indv_id, LAB_ID = wetlab_id, SAMPLE = sample_id, institution, 
                  RUN = batch, TISSUE = source_of_RNA, STATUS = affected_status, family_id, related_proband, 
                  age, sex, origin, disease_category, HPO_terms_ID, resolved_case, causal_gene, 
                  RIN, `260_280`, `260_230`, seq_machine)
  #print(head(as.data.frame(metadata_selected)[1:4,1:3]))
  # reformat for LIMMA/DREAM reteaped measures diff exp analysis
  metadata_reformat = metadata_selected %>% 
    mutate(build1 = paste0(SAMPLE, "_",builds[1]),
          build2 = paste0(SAMPLE, "_",builds[2])) %>% 
    gather(build, study_id, build1, build2) %>% 
    arrange(RUN, study_id)

  metadata_reformat_filt=metadata_reformat%>% filter(TISSUE == !!tissue) 
  metadata = column_to_rownames(metadata_reformat_filt, "study_id") 
   return(metadata)

}

load_count_data <- function(build, infile){
  #build_id = ifelse(grepl("chm", build), paste0(build, "ensembl"), build)
  #infile = paste0(input_dir, build, "/eOutliers/", tissue, ".", build_id, ".dedupOptical_minMQ255_rsem.genes.results")
  cat("\n\n")
  cat(paste(Sys.time(), build, "\n"))
  cat(paste0("  Reading file: ", infile, "\n"))
  
  df = fread(infile) %>% 
    rename(gene_id = V1) %>%
    # remove liftover indicator from gene_id
    mutate(gene_id = str_remove(gene_id, "_\\d+")) %>%
    # remove gene version from gene
    mutate(gene = str_remove(gene_id, "\\.\\d+"))

  # setnames(df,c("sample","gene_id","transcript_id","expected_count")) #set names is intentional: for size of df
  
  cat(paste("  loaded\n"))
  
  # df = df %>% 
  #   # create clean gene name that drops version but retains PAR annotations
  #   mutate(gene_id = str_remove(gene_id, "_\\d+")) %>%
  #   mutate(gene = str_remove(gene_id, "\\.\\d+")) %>% 
  #   dplyr::select(-gene_id, -transcript_id)
  
  # Filter out ensembl IDs that appear more than once per sample
  # (this is the result of mapping chm13 genes back to ensembl IDs when some  
  #  chm13 genes map to same ensembl ID, we exclude these from diff exp analysis)
  # OLD
  # get_dups = df %>% dplyr::filter(sample == df$sample[1]) # filters to first sample
  # dup_genes = get_dups[which(duplicated(get_dups$gene_id)), "gene_id"] %>% pull(gene_id) %>% unique() # for first sample, get list of repeated genes
  # cat(paste0("Removing duplicated genes: ", length(dup_genes)))
  # df = df %>% dplyr::filter(!(gene_id %in% dup_genes)) # drop repeated genes

  # for gene matrix
  # gene_vers_duplicates = as.data.frame(table(df$gene_id)) %>% filter(Freq > 1) %>% arrange(desc(Freq)) %>% rename(gene_id=Var1)
  # write_tsv(gene_vers_duplicates, paste0(output_dir1, tissue, ".", build , ".dedupOptical_minMQ255_rsem.genes.repeat_gene_versions.txt"))
  cat(paste(Sys.time(), "Checking for duplicate genes\n"))
  gene_duplicates = as.data.frame(table(df$gene)) %>% filter(Freq > 1) %>% arrange(desc(Freq)) %>% rename(gene=Var1)
  gene_duplicates_outfile = paste0(output_dir1, tissue, ".", build , ".dedupOptical_minMQ255_rsem.genes.repeat_genes.txt")
  write_tsv(gene_duplicates, gene_duplicates_outfile)
  cat(paste("  > Duplicate genes for", build, ":", nrow(gene_duplicates), "\n"))
  cat(paste("  Duplicate gene list saved for", build, ":", gene_duplicates_outfile, "\n"))

  # drop genes with multiple counts
  df = df %>% 
    filter(gene %notin% gene_duplicates$gene) %>%
    # dplyr::select(-gene_id) %>%
    # keeping gene_id instead of gene so that we can get the list of genes with different versions
    dplyr::select(gene, gene_id, everything())

  # NO LONGER NEEDED
  # df_wide = df %>% 
  #   pivot_wider(names_from = sample, values_from = expected_count) %>%
  #   # create clean gene name that drops version but retains PAR annotations
  #   mutate(gene_id = str_remove(gene_id, "_\\d+")) %>%
  #   mutate(gene = str_remove(gene_id, "\\.\\d+")) %>%
  #   dplyr::select(-transcript_id, -gene_id) %>%
  #   # dplyr::select(-transcript_id) %>%
  #   distinct()
  
  # append build to all columns (bc its fast), then remove from first column (gene_id)
  # so that only sampleIDs and gene versions have the build suffix
  setnames( df, paste0(colnames(df[1,]), "_", build) )
  setnames( df, paste0("gene_", build), "gene")
  
  cat(paste(Sys.time(), "Count matrix loaded for", build))

  return(df)
  
}



get_count_mtrx <- function(builds=c(build1, build2), infiles = c(build1_infile, build2_infile)){
  build1_counts = load_count_data(build = builds[1], infile = infiles[1])
  build2_counts = load_count_data(build = builds[2], infile = infiles[2])
  
  # cat("\n\n")
  # print(build1)
  # build1_counts[1:5,1:5]
  # cat("\n\n")
  # print(build2)
  # build2_counts[1:5,1:5]
  
  # get discordant gene versions
  cat("\n\n")
  cat(paste(Sys.time(), "Checking gene versions\n"))
  diff_gene_versions = full_join(build1_counts %>% dplyr::select(starts_with("gene")),
                                 build2_counts %>% dplyr::select(starts_with("gene"))) %>%
                                 filter(.[[2]] != .[[3]]) %>%
                                 distinct()
  diff_gene_versions_outfile = paste0(output_dir1, tissue, ".", build1, "_vs_", build2 , ".dedupOptical_minMQ255_rsem.genes.diff_gene_versions.txt")
  write_tsv(diff_gene_versions, diff_gene_versions_outfile)
  cat(paste("  > Number of genes with discordant gene versions:", nrow(diff_gene_versions_outfile), "\n"))
  cat(paste("  Discordant gene version list saved for", build1, "vs", build2, ":", diff_gene_versions_outfile, "\n"))

  # merge into counts matrix
  cat("\n\n")
  cat(paste(Sys.time(), "Merging count matrices\n"))

  ## check that samples are the same
  if( length(colnames(build1_counts)) != length(colnames(build1_counts)) ){
    stop(paste0("Column counts are not equal:\n  ", 
                build1, ": ", length(colnames(build1_counts)), 
                build2, ": ", length(colnames(build1_counts))))
  } else if( !setequal(colnames(build1_counts), colnames(build1_counts)) ){
    stop(paste0("The following samples are not paired:\n  ", 
                build1, ": ", colnames(build1_counts)[colnames(build1_counts) %notin% colnames(build2_counts)], 
                build2, ": ", colnames(build2_counts)[colnames(build2_counts) %notin% colnames(build1_counts)]))
  }

  count_matrix = inner_join(build1_counts %>% dplyr::select(-starts_with("gene_id")), 
                            build2_counts %>% dplyr::select(-starts_with("gene_id")), by="gene")

  count_matrix_clean = count_matrix %>% 
                          column_to_rownames(var = "gene") %>% 
                          dplyr::select(-starts_with("gene"))

  cat(paste(Sys.time(), "Count matrices merged\n"))

  return(count_matrix_clean)
  }

# Load count matrix
countMatrix = get_count_mtrx()

# save count matrix
write_tsv(countMatrix, opt$counts_outfile)
cat(paste(Sys.time(), "Saving merged count matrix:", opt$counts_outfile, "\n"))

# examine count matrix
n_samples = ncol(countMatrix)/2
cat(paste("  > n genes:", nrow(countMatrix), "\n"))
cat(paste("  > n samples:", n_samples, "\n"))
# dim(countMatrix) 
# countMatrix[1:5,1:5]
# countMatrix[1:5,(ncol(countMatrix)-4):ncol(countMatrix)]


# Load metadata
cat("\n\n")
cat(paste(Sys.time(), "Loading metadata\n"))

metadata1 = get_metadata(my_md, builds = builds, tissue)
# fwrite(metadata1, file="/oak/stanford/groups/smontgom/raungar/UDN/metadata.tmp", row.names=T)
# print(head(metadata1[1:5,1:5]))

# reorder metadata to match count matrix sample order
cat(paste0("  Reordering metadata to match count matrix sample order\n"))
metadata <- metadata1 %>% mutate(tmp_sample=rownames(metadata1)) %>%
                  dplyr::filter(tmp_sample%in%colnames(countMatrix[1,])) %>%
                    arrange(factor(tmp_sample,levels=colnames(countMatrix[1,])))

rownames(metadata)=metadata[,"tmp_sample"]
metadata$tmp_sample<-NULL
dim(metadata)
metadata[1:7,1:5]


## preprocessing ----

# filter genes by number of counts
cat("\n\n")
cat(paste(Sys.time(), "Filtering count data\n"))

# genes must have > 0  raw counts in at least 30% of the cohort
# genes must have > 0.1 cpm in at least 30% of the cohort
min_counts = 0
min_cpm = 0.1
min_prop = 0.3
cat(paste0("  Expression filters:\n", 
              "    greater than ", min_counts, " read counts and\n", 
              "    greater than or equal to ", min_cpm, " cpm (counts per million)\n", 
              "    in at least ", min_prop*100, "% of samples in each build (", 
              ceiling(min_prop*n_samples), " per build)\n\n"))

countMatrix=countMatrix[ sapply(countMatrix, is.numeric)] 
countMatrix=na.omit(countMatrix)
# print(countMatrix[1:7,1:3])

# requires a gene be expressed in 30% of samples across *either* build
# isexpr = rowSums(cpm(countMatrix)>0.1 & countMatrix>0,na.rm=T) >= 0.30*n_samples

# requires a gene be expressed in 30% of samples in *both* builds
counts_build1 = countMatrix[,1:(n_samples)]
counts_build2 = countMatrix[,(n_samples+1):(n_samples*2)]
# isexpr = rowSums(cpm(counts_build1 ) >= min_cpm & counts_build1 > min_counts &
#                     cpm(counts_build2 ) >= min_cpm & counts_build2 > min_counts,na.rm=T) >= min_prop*n_samples

isexpr1 = rowSums(cpm(counts_build1) >= min_cpm & counts_build1 > min_counts) >= min_prop*n_samples
isexpr2 = rowSums(cpm(counts_build2) >= min_cpm & counts_build2 > min_counts) >= min_prop*n_samples
isexp = cbind(as.data.frame(isexpr1), as.data.frame(isexpr2))
expressed = isexp %>% filter(isexpr1 & isexpr2) %>% rownames(.)
cat(paste0("  > Number of genes passed expression filter:", length(expressed), "\n"))

# save list of filtered genes
insufficient_exp = isexp %>% filter(isexpr1 == FALSE | isexpr2 == FALSE) %>% rownames_to_column(var = "gene")
setnames(insufficient_exp, c("gene", paste0("exp_", build1), paste0("exp_", build2)))
insufficient_exp_outfile = paste0(output_dir1, tissue, ".", build1, "_vs_", build2 , 
                                  ".dedupOptical_minMQ255_rsem.genes.not_enough_exp.", 
                                  min_counts,"counts_", min_cpm, "cpm_", min_prop*100, "pct_samples.txt")
cat(paste0("  > Number of genes failed expression filter: ", nrow(insufficient_exp), "\n"))
cat(paste0("      > in ", build1, ": ", insufficient_exp %>% filter(.[[2]] == FALSE & .[[3]] == TRUE) %>% nrow(), "\n"))
cat(paste0("      > in ", build2, ": ", insufficient_exp %>% filter(.[[2]] == TRUE & .[[3]] == FALSE) %>% nrow(), "\n"))
cat(paste0("      > in both: ", insufficient_exp %>% filter(.[[2]] == FALSE & .[[3]] == FALSE) %>% nrow(), "\n"))

cat(paste("  Saving list of genes with insufficient expression:", insufficient_exp_outfile, "\n"))
write_tsv(insufficient_exp, insufficient_exp_outfile)

cat("\n")
print(countMatrix[1:7,1:3])

rm(counts_build1, counts_build2, insufficient_exp)


## LIMMA/DREAM ----
max_pval = 0.01
min_logFC = 1

cat("\n\n")
cat(paste(Sys.time(), "Preparing for LIMMA-DREAM\n"))
# Normalization factors
# cat("\n")
cat(paste(Sys.time(), "Calculating normalization factors\n"))
geneExpr = DGEList( countMatrix[expressed,] )
geneExpr = calcNormFactors( geneExpr )


# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(1, progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ build + (1|SAMPLE) 

# print(rownames(metadata[,1]))
# print(dim(metadata))
# print(geneExpr[1:5,1:5])
# print(colnames(geneExpr[1,]))
# print(dim(geneExpr))

# estimate weights using linear mixed model of dream
# cat("\n")
cat(paste(Sys.time(), "Running voom with DREAM weights\n"))
vobjDream = voomWithDreamWeights( geneExpr, form, metadata,BPPARAM=param )
cat(paste(Sys.time(), "Voom complete\n"))
rm(metadata1, expressed)

# Fit the dream model on each gene
# cat("\n")
cat(paste(Sys.time(), "Fitting the differential expression model\n"))
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata, BPPARAM=param )
# fitmm = dream( geneExpr, form, metadata )
#fitmm_eb = eBayes(fitmm)
cat(paste(Sys.time(), "Differential expression calling complete\n"))

# Examine design matrix
# cat("\n")
cat(paste(Sys.time(), "Differential expression results\n"))
cat(paste0("BH-adjusted p-val <= ", max_pval, "\n",
             "absolute(log2FC) >= ", min_logFC, "\n"))

#outdir="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/TPMDifferences/limma/"

# coefs = coef(fitmm)
# coef_outfile=paste0(outdir, tissue,".dedupOptical_minMQ255.diff_expression_", build1, "_vs_", build2, ".limma_voom_dream.coefs.txt.gz")
# try(write_tsv(results,coef_outfile), silent=T)

# Get results of hypothesis test on coefficients of interest
topTable( fitmm, coef=paste0("build","build2"), number=5 ) # default adjustment is BH
results = topTable( fitmm, coef=paste0("build","build2"), n=Inf) %>% 
  as.data.frame() %>% 
  rownames_to_column("gene")
# dim(results)
# head(results)

n_sig = nrow(results %>% filter(adj.P.Val <= max_pval) %>% filter(abs(logFC) >= min_logFC))
cat(paste0("  > Total genes tested: ", nrow(results), "\n"))
cat(paste0("  > Significant genes: ", n_sig, "\n"))
cat(paste0("  > % Total that are significant: ", round(n_sig/nrow(results)*100,2), "\n"))

# results_outfile = paste0(outdir, tissue,".dedupOptical_minMQ255.diff_expression_", build1, "_vs_", build2, ".limma_voom_dream.all.txt.gz")
cat(paste("  Saving all results to:", results_outfile, "\n"))
write_tsv(results, results_outfile)

#  Sig results
cat(paste("  Saving significant results to:", opt$sig_outfile, "\n"))
sig = results %>% filter(adj.P.Val <= max_pval) %>% filter(abs(logFC) >= min_logFC)
write_tsv(sig, opt$sig_outfile)

# Fitmm
head(fitmm$design, 3)
cat(paste("  Saving fit model to:", opt$fitmm_outfile, "\n"))
write.fit(fitmm, NULL, opt$fitmm_outfile, adjust="BH", method ="separate")


## Make Contrasts ----

# head(coefs)
# contr <- makeContrasts(build - build2, levels = colnames(coef(fitmm)))
# tmp <- contrasts.fit(fitmm, contr)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)


