
conda activate phen2gene
mkdir -p ./gene_lists/

for i in ./hpo_lists/*; do 
    file_name=$(basename $i)
    sample=${file_name%%.*}
    python3 phen2gene.py -f $i -out ./tmp 
    head -251 ./tmp/output_file.associated_gene_list > ./gene_lists/${sample}.phen2gene.gene_list.txt #get top 250 prioritized genes + 1 header line
    rm -rf ./tmp

    cat ../gene_lists/${sample}.phen2gene.gene_list.txt >> ../gene_lists/merged_gene_list.txt
done
