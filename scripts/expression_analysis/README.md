# Expression Outlier Calling
(note: Snakefile is an old script for reference, will be used to build upon. The main Snakefile to use is SnakefileSingleSamplesheet)      

## SnakefileSingleSamplesheet   
get_rsem_by_tiss
 * input: all RSEM files     
 * params: tissue  
 * delimeter to extract on   
 * output: per tissue RSEM combined files     
 * purpose: combine all files from a single tissue      

get_eoutlier_gene_or_isoform    
 * input: custom outlier calling script, metadata, tissue combined file    
 * output: z-scores as well as intermediate files for batch checking     
 * purpose: call outliers and get intermediate files    

split_by_sample   
 * input: per tissue file    
 * output: per individual files   
 * purpose: get z scores at a per-individual level for all sampels  


## Other files
manually_test_batch.Rmd    
 * this script is manual, used to identify batch correction and look at general QC    
get_med_tpm.R   
 * this is used for later on, GCs wanted to know median TPM from GTEx for a gene by tissue   
split_expr_by_tissues.py     
 * this is used to extract the GTEx reference data   
TPM_to_[].R   
 * all of these need their own script due to metadata differences   
