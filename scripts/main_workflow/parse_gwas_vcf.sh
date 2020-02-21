# load bcftools on BC4
eval "module load apps/bcftools-1.9-74/1.9-74"

data="/XremovedX/"
ids_list="$1"  # input column file of trait IDs

while read i; do 

# extract required fields from VCF
cmd="bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP]\n' ${data}/${i}/${i}.vcf.gz -o tmp/${i}_subset.tsv"
echo $cmd
eval $cmd

# keep the required cols only 
cmd="less tmp/${i}_subset.tsv | cut -f3-8 > tmp/${i}_subset_slim.tsv"
echo $cmd
eval $cmd


# add column names 
cmd="(echo -e 'rsid\tallele_0\tallele_1\t${i}_b\t${i}_se\t${i}_pval'; cat tmp/${i}_subset_slim.tsv) > S_XY2/${i}.tsv"
echo $cmd
eval $cmd


done < $ids_list
