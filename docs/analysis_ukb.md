# Case study 1: UK Biobank, setting up metaCCA analysis

Current location: project + `working/data/`

Plink (on epi-franklin): `$plink = /data/ny19205/software/plink/plink`

## Selecting traits (IEU-GWAS db)

See `scripts/select_traits`

To select traits I used the reference table `study-table-13-02-19_UKBonly.tsv` from UKB data (`biobank_traits_parser.R`).


UKB traits:
	
	* BMI: 19953 
	* Waist cirmunverence: 9405
	* Whole body fat mass: 19393
	* Whole body fat-free mass: 13354
	* Body fat percentage: 8909
	* Systolic blood pressure (automated): 20175
	* Diastolic blood pressure (automated): 7992

UK biobank SNP list (the same for all traits):  9,837,128

## Preparing the input data

The data for selected IEU UKB studies is available in VCF format here: `/mnt/storage/private/mrcieu/research/scratch/IGD/data/public/` (BC4)

UKB IEU data was parsed using bash script `parse_gwas_vcf.sh`, which can be used on other datasets available in VCF. The script was run on BC4 in `/mnt/storage/scratch/ny19205/`, then clean tsvs were moved to project space to `data/S_XY_matrices/ukb` and `data/S_XY_matrices/ukb`. 

#### Reference data

1. Create a list of SNPs in UKB

	```bash
	cd snps_list/
	
	less UKB-b-13354_subset.tsv | cut -f3,4,5 | sed 's/\t/_/g' | sort | uniq > ukb_common_snps.txt
	
	# save as range	
	cat ukb_common_snps.txt | sed 's/_/\t/g'| awk '{print $1"\t"$2"\t"$2"\t"$3}' > ukb_common_snps_chrpos.txt
	wc -l ukb_common_snps_chrpos.txt
		9837128
		
	cat ukb_common_snps.txt | cut -f3 -d"_"| sort| uniq > ukb_common_snps_rsid.txt
	wc -l ukb_common_snps_rsid.txt
		9834852	
	```


2.  Subset reference data to the GWAS list of SNPs
	
	```bash
	cd ../genotype_matrix_2/
	plink --bfile ../1000GPdata/data_maf0.01_rs_ref --extract ../snp_lists/ukb_common_snps_rsid.txt --exclude ../snp_lists/sticky_snps.txt --out data_overlap --make-bed --keep-allele-order --snps-only
	wc -l data_overlap.bim
		7,565,603
	```
	what this does:

	* subset only to position in UKB
	* drop sticky SNPs (they create NaN in LD calculation)
	* drop SNPs that have probelms with alleles
	* restrict flipping of alleles
	* keep only SNPs (not indels)

## Creating genotype-phenotype matrices (S_XY) for studies

To join the within-study trait data used script `1_prepapre_data_XY.Rmd` which saves the data in `geno-pheno_matrix_UKB.RData`.


## Creating genotype matrix (S_XX) from reference data


#### Pruning data
```bash
plink --bfile data_overlap  --indep-pairwise 50 5 0.01 --out data_overlap
plink --bfile data_overlap --extract data_overlap.prune.in --make-bed --out data_pruned  --keep-allele-order
wc -l data_pruned.bim
	139,735
```

#### Annotate SNPs with gene names
```bash
cat <(echo -e "CHR\\tSNP\\tBP\\tREF\\tALT") <(cut -f1,2,4,5,6 data_pruned.bim) > data_pruned_snps.txt
plink --bfile --annotate data_pruned_snps.txt ranges=../snp_lists/glist-hg19 --out data_pruned
	60418 out of 139735 annotated
```

#### Subset data to the annotated genes 
```bash
# first select SNPs that got annotated with gene names
less data_pruned.annot| tr ' ' \\t |  grep -v "\." > annotated_genes.txt
	60361 # SNPs
	
# output SNP list
less annotated_genes.txt | grep -v "^CHR"| cut -f2 > annotated_SNPs_list.txt	
	
# count SNPs per gene
less annotated_genes.txt | cut -f6 | sort | uniq -c | sort -k1n > annotated_SNPs_per_gene.txt
	13531 # genes
```
13,531 genes -- 60,361 SNPs -- 1-227? SNPs per gene -- on average 4.47? SNPs per gene

#### Subset data to SNPs to include, and split into chromosomes

```bash
plink --bfile data_overlap  --extract annotated_SNPs_list.txt  --out data_overlap_subset --keep-allele-order --make-bed
```

```bash
mkdir data_overlap_subset_by_chr
for chr in {1..22}; do
$plink --bfile data_overlap_subset --chr $chr --out data_overlap_subset_by_chr/data_overlap_chr${chr} --keep-allele-order --make-bed; 
done
 
# view SNPs by chr
wc -l data_overlap_subset_by_chr/*bim | sort -k1n
```


#### Create LD matrix from SNPs that belong to genes
Creates an *LD matrix of r or r2 values*  from 502 European samples from 1000 Genomes phase 3 data.

```bash
mkdir ld_matrix_by_chr2/
for chr in {1..22}; do 
$plink --bfile data_overlap_subset_by_chr/data_overlap_chr${chr} --r2 square --out ld_matrix_by_chr/ld_matrix_chr${chr} --keep-allele-order;
 done
 
mkdir ld_matrix_by_chr/
for chr in {1..22}; do 
$plink --bfile data_overlap_subset_by_chr/data_overlap_chr${chr} --r square --out ld_matrix_by_chr/ld_matrix_chr${chr} --keep-allele-order;
 done 
```

The 22 `ld_matrix_chrN.ld` files (10-379 MB) are the LD matrces that contain squared Pearson corrlation coefficients between the selected SNPs in the chromosome. Each matrix will be used as input to the metaCCA method as genotype-genotype matrix. 

#### Tidying up LD matrix

Add row and column names (SNPs) to each chr LD matrix in R, to make it suitable to be used as S_XX in metaCCA. 

See script: `1_prepapre_data_XX_per_chr.Rmd`

The output is `S_XX_matrices/LDmatrix_chrN.Rdata`  files that can be read in by the main metaCCA script.




## Running metaCCA

Currently simple testing is done in `2_analyse_data_testing.Rmd `. 

Both

- Univariate SNPs - multiple phenotypes analysis 
- Multivariate SNPs - multiple phenotypes analysis 

are scripted in `3_run_metaCCA_analysis.R` to be run in parallel process on BC3.

Submission script on BC3 (or here) in `/newhome/ny19205/metaCCA_big_jobs`: `submit_metaCCA_jobs.sh` 

Submit one chr job to one node (16 cores):

```qsub -l walltime=05:00:00 submit_metaCCA_jobs.sh -F "22"```

## Post-metaCCA results processing

### Multivariate SNPs (gene-based analysis)

The output is 22 files `metaCCA_multisnp_genes_chrN.tsv`

Merge then and tidy up like this:

```bash
cat metaCCA_multisnp_genes_chr* | gsed  -z 's/,\n/,/g' | gsed  -z 's/\n"/"/g'  |  sort -k2 -r | uniq > multivar_metaCCA_merged_results.tsv
``` 

Not much intresting after that at the moment.

### Univariate SNPs (SNPs-based analysis)

`4_review_results_gwascat.Rmd` - merge results with GWAS catalog and annotate by groups

In the middle of this R script can get LD proxies using Python script `python_LDproxies/get_LD_proxies.py` to match more things in GWAS catalog

`5_visualise.Rmd` - boxplots, dotplots, barplots


### GWAS catalog annotation plots

<img src="../plots/UKB_9M/boxplot.png" alt="diagram" width="500"/><br>
<img src="../plots/UKB_9M/dotplot.png" alt="diagram" width="500"/>