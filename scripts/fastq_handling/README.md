# Process new batch of RNA-seq data for rare disease project
This first trims adaptors using cutadapter, makes a reference file if not created for STAR, then runs STAR in a two-pass manner to map reads. Finally, reads are filtered for high enough mapping quality and PCR duplicates are removed using Picard. At last, RSEM is used to quantify the reads. This is updated for hg38.    
This needs to be slightly manually modified for each batch for the proper batch number and sample names


## Order
#### SnakefileSetupFASTQ   
This snakefile only needs to be run once, it preps STAR and RSEM files and is not sample dependent
#### SnakefilePreAlignment  
This snakefile is run on the FASTQ files post-bcl processing  
#### Snakefile   
This snakefile aligns and quantifies reads   



## SnakefileSetupFASTQ   
#### This file prepares the rsem reference and the STAR index. This is required input for running STAR and RSEM

### Required config file keys:
* `build` - name of reference genome build; STAR indices will be labeled with this for posterity
* `gtf_unzip` - path to uncompressed gtf (genome annotation)
* `fa_unzip` - path to uncompressed fasta file (genome sequence)
* `rsem_ref` - path to directory to store RSEM references in
* `rsem_dir` - path to directory to store RSEM reference creation log in
* `reference_files_dir` - path to directory to store STAR references in
* `star` - path to the STAR executable

### Rules
`rule make_rsem_reference`:
 * purpose: make a reference file for downstream RSEM
 * input: gtf, fasta
 * output: RSEM reference file (location and prefix specified in config)
 * *to_skip:* remove `config["rsem_dir"]+"/rsem_reference.log"` from `rule all`

`rule create_star_annot_75`:   
 * purpose: create STAR index annotation file for reads that are length 76    
 * input: fasta, gtf    
 * output: log file, index file for future use (location and prefix specified in config file) 
 * *to_skip:* remove `config["reference_files_dir"]+"/STAR2.7.4a_INDEX_"+build+"_ov75.log"` from `rule all`    

`rule create_star_annot_99`:    
 * purpose: create STAR index annotation file for reads that are greater than length 99   
 * input: fasta, gtf  
 * output: log file, index file for future use (location and prefix specified in config file)
 * *to_skip:* remove `config["reference_files_dir"]+"/STAR2.7.4a_INDEX_"+build+"_ov99.log"` from `rule all` 

## SnakefilePreAlignment
#### The purpose of this file is to prep the FASTQs. This is independent of alignment, so has been put in its own snakefile.

> Note: Still working on making the demultiplexing step optional for pipelines starting from untrimmed fastq files, and creating an environment with all required software.

### Required config file keys:
* `sample_file` - path to sample table
* `bcl2fastq_dir` - path to uncompressed gtf (genome annotation)
* `fastq_dir` - directory for fastq files; structure: `${fastq_dir}`/`FASTQ_RUN${run}`/`sample1.fastq`
* `qc_dir` - directory to save quality control reports to
* `fastqc` - path to fastqc executable

### Required sample table columns (case-insensitive):
* `SAMPLE` - [required] unique identifyers for each sample
* `BCL_DIR` - [optional] path to the BCL data that needs to be demultiplexed. This directory must contain a file called Samplesheet.csv that contains the demultiplexing information (including sample names and indices). This file is typically generated by the sequencing machine software.
* `RUN` - [optional] Identifier for sequencing run; If RUN is not specified, all outputs will be written to `config["fastq_dir"]/FASTQ_RUN`
* `READ_LEN` - [required] specify sequencing read length per sample
* `INDEX1` - [required] i5 index sequence 1 (can be obtained from Samplesheet.csv in the raw BCL sequencing directory)
* `INDEX2` - [required] i7 index sequence 2 (can be obtained from Samplesheet.csv in the raw BCL sequencing directory)

### Rules
`rule trim`:   
 * purpose: trims adaptors with cutadaptor    
 * input: SampleSheet.csv and read1.fastq.gz and read2.fastq.gz   
 * output: log file used downstream, read1.trimmed.fastq.gz, read2.trimmed.fastq.gz  

`rule zip`:   
 * purpose: gzip the trimmed files
 * input: trimmed fastqs
 * output: gzipped trimmed fastqs    

`rule run_fastqc`:
 * purpose: get qc statistics on these trimmed files  
 * input: fastqc tool, sample of fastq
 * output: fastq html report     

## Snakefile

> Note: Still working on creating an environment with all required software.

### Required config file keys:
* `build` - name of reference genome build (should match the label used when creating the STAR indices)
* `sample_file` - path to sample table
* `fastq_dir` - directory for fastq files; structure: `${fastq_dir}`/`FASTQ_RUN${run}`/`sample1.fastq`
* `star_dir` - directory for STAR output to be stored
* `bam_dir` - directory for filtered bams
* `rsem_dir` - directory for RSEM quantification results to be stored
* `qc_dir` - directory to save quality control reports to
* `fastqc` - path to fastqc executable
* `rsem` - path to the rsem executable
* `samtools` - path to samtools excutable
* `picard` - path to picard executable
* `star` - path to the STAR aligner executable
* `star_annot_ov75` - path to the STAR annotation files for samples with reads <100bp
* `star_annot_ov99` - path to the STAR annotation files for samples with reads >99bp

### Required sample table columns (case-insensitive):
* `SAMPLE` - [required] unique identifyers for each sample
* `RUN` - [optional] Identifier for sequencing run; If RUN is not specified, all outputs will be written to `config["fastq_dir"]/FASTQ_RUN`
* `READ_LEN` - [required] specify sequencing read length per sample
* `VCF_FOR_ASE` - [optional] VCF containing biallelic variants for to include in a WASP implementation of STAR to correct for allelic mapping bias; if no VCF is provided, STAR will run without WASP

### Rules
`rule run_star`:
 * purpose: run STAR
 * input: star annotation, fastq out directory, r1 fastq, r2 fastq, log file
 * output: completed star file
 
`rule mapping_filter`:
 * purpose: filter for mapping quality greater than 30 and then sort
 * input: star output file
 * output sorted bam file

`rule rm_pcr_duplicates`:
 * purpose: to remove PCR duplicates using picard
 * input: sorted bam
 * output: dedup bam, metrics log file


`rule convert_bam_for_rsem`:
 * purpose: to allow rsem to modify bam for proper quantification
 * input: bam file from rule rm_pcr_duplicates
 * output: dedup_rsem.bam file

`rule run_rsem`:
 * purpose: actually run RSEM to quantify
 * input: sorted and dedupped bam file  
 * output: log file, and the rsem file    

`rule run_multiqc`:   
 * purpose: for qc analyses and understanding downstream   
 * input: multiqc tool, fastqc directory, star direcotry, bam directory, 
 * params: prefix for multiqc report and outdir      
 * output: fake log file since output is not predictable       

### SnakefileSINGLEEND
For DGN, a new Snakefile was created to specifically deal with the fact it is single-ended. All rules, except the first one, shoudl be the same as previously, just with specific parameters changed for processing single-ended reads          
Therefore, the only rule explained will be the new one    