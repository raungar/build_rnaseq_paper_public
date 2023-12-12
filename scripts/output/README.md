## Creating Output Files for GCs    
The point of this Snakefile is to create the output for GCs.    

combine_everything    
 * input: everything you would like to input, NAs where they do not exist (can look through and identify that), Rscript    
 * output: per sample combined file   
 * purpose: this is the meat of the snakefile. this uses a custom script to combine all fields and format them for GCs    

convert_for_gcs  
 * input: combined file      
 * output: renamed combined file    
 * purpose: this just renames the file as per request by GCs    
