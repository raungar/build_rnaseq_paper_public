# Perform splicing outlier analysis
This uses LeafcutterMD and Regtools to perform the splicing outlier analysis


index_bam_for_regtools  
 * input: bam file aligned to genome   
 * output: bam.bai indexed bamfile   
 * purpose: bam needs to be indexed for regtools junction calling   

bam_to_junc  
 * input: bam file aligned to genome,bam2junc.sh from leafcutter  
 * output: junction file   
 * purpose: extract junctions (using leafcutter original script)

bam_to_junc_regtools   
 * input: indexed bam and bam aligned to genome, regtools   
 * output: junction file   
 * purpose: extract exon exon juctions   

annotate_juncs     
 * input: junction file, reference fasta, reference gtf, regtools   
 * output: annotated junction file via regtools    
 * purpose: get more information about these exon-exon junctions    

get_all_samples   
 * input: all junc files, metadata file   
 * output: combined junction files across a tissue    
 * purpose: for each tissue, combine all the junction files for that tissue to input into leafcutterMD   

intron_clustering    
 * input: tissue-specifc junctions, leafcutter script    
 * params: min_read_support, max_intron_len, output directory, prefix for these files   
 * output: dir/prefix_perind.counts.gz and dir/prefix_perind_numers.counts.gz   
 * purpose: leafcutter step to cluster introns together, see leafcutter documentation for more details    

outlier_calling   
 * input: script from from leafcutterMD, numers.counts.gz file   
 * params: prefix for output    
 * output: pVals file     
 * purpose: calculate pvalues from these intron clusters, see leafcutterMD documentation for more details  

p_to_z    
 * input: script from  leafcutterMD, pvals from previous rule    
 * output: zscore file     
 * purpose: multiple test correct these pvals, and then convert to z-scores    

split_by_sample   
 * input: zscore file   
 * params: outdir prefix   
 * output: file per individual   
 * purpose: split file by individual for individual-level reports

sort_split   
 * input: split file by individual, bedtools      
 * output: sorted file from above   
 * purpose: sort for annotation step    

zannotate_juncs
 * input: annotate_juncs, sorted file, bedtools   
 * output: annotated individual file    
 * purpose: annotate z scores with junc level information

annotate_juncs_genes (i want to change this to use regtools)   
 * input: gtf reference, zannotate_juncs output, bedtools   
 * output: annotated individual file    
 * purpose: actually include the nearby genes

uniq_feat_per_anno (will change for regtools)    
 * input: annotated by individual   
 * output: annotated by individual reduced  
 * purpose: deals with multiple annotations.   

collapse_top_feature     
 * input: custom script, previous reduced annotation file    
 * output: collapsed file    
 * purpose: collapsed at the transcript or gene level   

combine_features    
 * input: all collapsed sample files    
 * output: combined sample files    
 * purpose:  combine all these files together for downstream analyses     
