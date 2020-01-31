
library(metaCCA)
library(plyr)
library(dplyr)
library(readr)
library(tibble)
library(purrr)
library(parallel)
library(MASS)
library(pryr)


memory_comments = T

if ( memory_comments == T ) print(paste0("Starting memory: ", pryr::mem_used()))

if(!exists("args") || class(args)!="character"){
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args)==0) {
  
  cat("MetaCCA analysis script
      
      Usage: Rscript <this_script>.R <numCores> <analysis_type> <chrs> [<data_path>]
      * <numCores> : number of cores to use to run n metaCCA analyses in paralell
      * <analysis_type> : 'snps' - univariate SNPs analysis; 
                          'genes' - multivariate SNPs per gene analysis (genes that contain 2+ SNPs);
      * <chrs> : chromosome to process; the parameter accepts individual number 1-22 or 'all' to do all chromosomes; if analysis snps specify 'all'
      * <data_path>: data directory; if not supplied, will assumed that project folder on epi-franklin is to be used
      ")
  stop("Must provide at least 3 arguments")
  
} else if ( length(args) == 3 ){
  
  numCores <- args[1] 
  analysis_type <- args[2] 
  chrs <- args[3]
  data_path <- "/projects/MRC-IEU/research/projects/ieu2/p4/041/working/data/" # this works on epif

  # data_path<-"/newhome/ny19205/metaCCA_big_jobs/data/" BC3
  
} else{
  numCores <- args[1]
  analysis_type <- args[2]
  chrs <- args[3]
  data_path <- args[4]
}

print(paste("Set cores:", numCores))

if(!analysis_type %in% c("genes", "snps")){
  stop("Not acceptable analysis type")   
} else{
  print(paste("Analysis type:", analysis_type))
}

if (!chrs  %in% c(1:22) & chrs != "all"){
  stop("You have incorrectly specified chromosomes to analyse; script accepts either individual chr number or 'all' ")
} else{
  if (chrs =="all"){chrs <- rev(1:22)} 
  print(paste("Selected chromosome(s):", chrs))
}

if (!dir.exists(data_path)){
  stop("Specified data directory does not exist")
} else{
  print(paste("Data location:", data_path))
}

if ( memory_comments == T ) print(paste0("Memory after loading libraries and arguments: ", pryr::mem_used()))


print("Loading genotype data...")
load(paste0(data_path, "S_XY_matrices/geno-pheno_matrix_UKB.RData")) #s_xy_ukb
N = 436419

if ( memory_comments == T ) print(paste0("Memory after loading S_XY: ", pryr::mem_used()))


print("Loading and parsing gene annotations to create gene dictionary... ")
# list of SNPs to exclude from the analysis
snps_to_drop<-read_tsv(paste0(data_path,"snp_lists/bad_allele_snps.txt"), col_names = F) %>% as_vector() %>% unname()

# SNPs in genes list
ref<-read_tsv(paste0(data_path,"genotype_matrix_2/annotated_genes.txt"), col_names=T) %>% filter(!SNP %in% snps_to_drop)

if ( memory_comments == T ) print(paste0("Memory after loading text files: ", pryr::mem_used()))





print("Creating S_YY matrix... ")
s_yy = estimateSyy( S_XY = s_xy )
if ( memory_comments == T ) print(paste0("Memory after estimating S_YY: ", pryr::mem_used()))


#### FUNCTIONS
apply_metaCCA_multi <- function(current_gene, gene_dict, s_xy, s_yy, s_xx, N) { 
  
  print(paste0("Running metaCCA on gene ", current_gene))
  print(paste0("It contains ", length(gene_dict[[current_gene]]), " SNPs: ")) 
  cat(gene_dict[[current_gene]])
  cat("\n\n")
  
  metaCCA_multi = metaCcaGp( nr_studies = 1, 
                             S_XY = list( s_xy ), 
                             std_info = c( 1 ),  #  standardised
                             S_YY = list( s_yy ),
                             N = c( N ) ,
                             analysis_type = 2,
                             SNP_id = gene_dict[[current_gene]],
                             S_XX = list( s_xx ) ) 
  
  metaCCA_multi_upd<-metaCCA_multi %>% 
    rownames_to_column("SNPs") %>% 
    mutate(SNPs= gsub('[()"c ]', '', SNPs)) %>% 
    #mutate(trait_weights= gsub('[()"c ]', '', trait_weights)) %>% 
    #mutate(snp_weights= gsub('[()"c ]', '', snp_weights)) %>% 
    mutate(gene = current_gene) %>% 
    dplyr::rename(log10pval=`-log10(p-val)`) %>% 
    mutate( pval = ifelse(is.infinite(log10pval), 0, 10^-log10pval)) %>%  # if logpval is Inf, convert it to 0 as pval
    mutate(gene = current_gene) %>% 
    mutate(snps_count = length(gene_dict[[current_gene]])) %>% 
    dplyr::select(c("gene", "r_1", "pval", "SNPs", "snps_count"))#, "trait_weights", "snp_weights"))
  
  result<-rbind(data.frame(), metaCCA_multi_upd)
  return(result)
}  


apply_metaCCA_snp <- function(current_snp, ref, s_xy, s_yy,  N) { 
  
  print(paste0("Running metaCCA on SNP ", current_snp))

  metaCCA_single = metaCcaGp( nr_studies = 1, 
                              S_XY = list( s_xy ), 
                              std_info = c( 1 ),  # standardised
                              S_YY = list( s_yy ),
                              N = c( N ) ,
                              analysis_type = 1,
                              SNP_id = current_snp) 
  
  metaCCA_single_upd<-metaCCA_single %>% 
    rownames_to_column("snp") %>% 
    mutate(snp= gsub('[()"c ]', '', snp)) %>% 
    #mutate(trait_weights= gsub('[()"c ]', '', trait_weights)) %>% 
    dplyr::rename(log10pval=`-log10(p-val)`) %>% 
    mutate(pval= 10^-log10pval) %>%
    left_join(ref[,c("snp", "ANNOT")], by= "snp") %>% 
    dplyr::rename(gene = ANNOT) %>% 
    dplyr::rename(SNP = snp) %>% 
    dplyr::select(c("SNP", "r_1", "pval", "log10pval", "gene" ))#,"trait_weights"))
    
  result<-rbind(data.frame(), metaCCA_single_upd)
  return(result)
}  
#####



# depending on the selected analysis type, do preparations
if ( analysis_type == "genes"){ 
  
  print(". . . Creating GENE:SNPs dictionary for multiSNPs(=gene)-multiple traits analysis ")
  
  # create gene:snps list
  # this is a very ugly way of doing it
  gene_dict<- setNames(vector("list", length(unique(ref$ANNOT))), unique(ref$ANNOT))
  ref$snp<-paste0( ref$SNP, "_", ref$REF, "_", ref$ALT)
  for (gene in unique(ref$ANNOT)){
    snps<- ref %>% filter(ANNOT == gene) %>% 
                    dplyr::select(snp) %>% 
                    as_vector() %>% unname()
    gene_dict[[gene]] <-snps
  }
  
  # create summary of SNPs per gene
  gene_dict_summary <- setNames(data.frame(matrix(ncol = 2, nrow = length(gene_dict)) ),
                                                      c("gene", "snp_count"))
  gene_dict_summary$gene <- as.character(gene_dict_summary$gene)
  gene_dict_summary$snp_count <- as.numeric(gene_dict_summary$snp_count)
  for (i in 1:length(gene_dict)){
    gene_dict_summary$gene[i] <- names(gene_dict)[i]
    gene_dict_summary$snp_count[i]<-length(gene_dict[[i]])
  }
  
  gene_dict_summary<-
    inner_join(gene_dict_summary, ref[, c("CHR", "ANNOT")], by=c("gene"="ANNOT")) %>% 
    unique() %>% 
    arrange(by_group=desc(snp_count))
  write_tsv(gene_dict_summary, paste0(data_path, "results/gene_snp_count_summary.tsv"))
  
  
  # select multip-SNP genes and single SNP genes
  #multi_snp_genes<-gene_dict_summary %>% filter(snp_count > 1) %>% 
  #                dplyr::select(gene) %>% as_vector() %>% unname()
  #
  #single_snp_genes<- gene_dict_summary %>%  filter(snp_count == 1) %>% 
  #                  dplyr::select(gene) %>% as_vector() %>% unname()
  #
  if ( memory_comments == T ) print(paste0("Memory after creating dicts: ", pryr::mem_used()))
  
  
} else if (analysis_type == "snps"){
  
  print(". . . Creating a list of SNPs for singleSNP-multi traits analysis ")
  # get the list of SNPs to analyse
  ref$snp<-paste0( ref$SNP, "_", ref$REF, "_", ref$ALT)
  single_snps<-as.vector(ref$snp)
}



## Run metaCCA for the selected analysis type

if (analysis_type == "genes"){
  
  cat("\n\n")
  print(paste("Running multiple SNP - multiple traits metaCCA analysis"))
  cat("\n\n")
  
  #chrs=c(13,7,2,16,9,6,4,3,1)
  for (chr in chrs ){
  
    print(paste0("Loading LD matrix for chr ", chr))
    load(paste0(data_path, "S_XX_matrices/LDmatrix_chr",chr ,".RData")) # s_xx
    
    if ( memory_comments == T ) print(paste0("Memory after loading S_XX ", pryr::mem_used()))
    
    # select gene that are only in current chr:
    gene_list <- gene_dict_summary %>% 
      filter(CHR== chr) %>% 
      filter(snp_count > 1) %>% 
      dplyr::select(gene) %>% 
      as_vector() %>% unname()
    
    metaCCA_multi_list <- mclapply(gene_list, apply_metaCCA_multi, mc.cores = numCores,
                                   gene_dict = gene_dict, 
                                   s_xy = s_xy, s_yy = s_yy, s_xx = s_xx,
                                   N = N)
    rm(s_xx) # remove from memory
    if ( memory_comments == T ) print(paste0("Memory after removing chr S_XX ", pryr::mem_used()))
    
    # convert list of vectors to df
    metaCCA_multi_df <- bind_rows(lapply(metaCCA_multi_list, as.data.frame.list))
    
    print("Saving multi SNP metaCCA data... ")
    write_tsv(metaCCA_multi_df, paste0(data_path, "results/metaCCA_multisnp_genes_chr",chr,".tsv"))
  }

  
  
} else if (analysis_type == "snps"){
  
  cat("\n\n")
  print(paste("Running single SNP - multiple traits metaCCA analysis"))
  cat("\n\n")

  metaCCA_single_snps <- mclapply(single_snps, apply_metaCCA_snp, mc.cores = numCores,
                                  ref = ref, 
                                  s_xy = s_xy, s_yy = s_yy, 
                                  N = N)
  
  # convert list of vectors to df
  metaCCA_single_snp_df <- bind_rows(lapply(metaCCA_single_snps, as.data.frame.list))
  print("Saving data single-SNP gene results... ")
  write_tsv(metaCCA_single_snp_df, paste0(data_path, "results/metaCCA_snps_all.tsv"))
  
  
}


  



