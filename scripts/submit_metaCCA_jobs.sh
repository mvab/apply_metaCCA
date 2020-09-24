#!/bin/bash

# input from command line

chr="$1"  # chromosome to process; add ' -F "chr#" ' at the end to submission script
cores=16
data_path="/newhome/XremovedX/metaCCA_big_jobs/data/"


# request resources:
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00

## not using: #PBS -q himem


# Define working directory
export WORK_DIR=$HOME/metaCCA_big_jobs

# load latest R
module add languages/R-3.5.1-ATLAS-gcc-6.1

# on compute node, change directory to 'submission directory':
cd $WORK_DIR


echo JOB ID: $PBS_JOBID
echo Working Directory: `pwd`
echo Start Time: `date`

# run your program, timing it:

# gene-based analysis
#time Rscript 3_run_metaCCA_analysis_IEU-NL.R  ${cores} genes $chr ${data_path}
# snp-based analysis
#time Rscript 3_run_metaCCA_analysis_IEU-NL.R  ${cores} snps all ${data_path}


# snp-based analysis NEW SCRIPT
#list1
#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/annotated_genes_list1_wo_dis_beta.txt ${data_path}

#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_GIANT_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/annotated_genes_list1_wo_dis_beta.txt ${data_path}



#list2
#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/list2_pval10e8_df.tsv ${data_path}

#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_GIANT_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/list2_pval10e8_df.tsv ${data_path}



#list3
#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/list3_pval10e4_df.tsv ${data_path}

#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_GIANT_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/list3_pval10e4_df.tsv ${data_path}



#list4
time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/list4_pval10e3_df.tsv ${data_path}

#time Rscript 3_run_metaCCA_analysis_upd.R  ${cores} snps all 339224 S_XY_matrices/geno-pheno_matrix_UKB_GIANT_2M_updGIANT.RData genotype_matrix_3/annotated_genes_subsetting/list4_pval10e3_df.tsv ${data_path}

echo End Time: `date`

