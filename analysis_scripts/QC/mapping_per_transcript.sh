#!/bin/bash
module load samtools 

dir_prefix="/oak/stanford/groups/smontgom/shared/UDN/Output"
suffix="Aligned.toTranscriptome.out.bam"
outfile_prefix="/oak/stanford/groups/smontgom/shared/UDN/Analysis/ReferenceComparisonPrimary/Mapping/mapping_"

for this_build in hg19 hg38 chm13
do
    outfile=$outfile_prefix"_"$this_build".txt"

    files=`ls $dir_prefix/$this_build/STAROutput/*$suffix`
    for file in $files 
    do
        ##get sample
        this_sample=LSJDflsdfjsf
        #parallelize
        samtools view $file | awk -F"\t" '{print $3, $(NF-2)}' | awk -F"[ |:]" -v build=$this_build -v sample=$this_sample '{print build"\t"sample"\t"$1"\t"$NF}' >> $outfile
        echo $file
        break
    done
    break
done