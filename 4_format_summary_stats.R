
setwd("/panfs/jay/groups/17/leej5/edwa0506/voting_gwas")


library(qqman)
library(vroom)
library(data.table)
library(tidyverse)

format_bolt <- function(bolt_file,  # Bolt summary statistics
                        pheno_file, # phenotype file with no missing samples
                        info_file,  # file with call rate and info
                        hwe_file,   # HWE p values per snp
                        output_file #
                        ) {
  
  columns <- c("SNPID",
               "CHR",
               "BP_b37",
               "EFFECT_ALLELE",
               "OTHER_ALLELE",
               "EAF",
               "BETA",
               "SE",
               "P",
               "N",
               "INFO",
               "HWE_PVAL",
               "CALLRATE")
  
  pheno <- vroom(pheno_file)
  sample_size <- nrow(pheno)
  
  HWE <- fread(hwe_file, stringsAsFactors = FALSE) 
  HWE <- HWE |> 
            select(SNP, P)|>
            rename(SNPID = 'SNP',
                   HWE_PVAL ='P')
  
  info <- vroom(info_file)
  info <- info |>
              select(rsid, AvgCall, Rsq) |>
              rename(SNPID = 'rsid',
                     INFO = 'Rsq',
                     CALLRATE = 'AvgCall')
  
  dat <- vroom(bolt_file)
  dat <- dat |> 
            rename(SNPID = 'SNP',
                   BP_b37 = 'BP',
                   EFFECT_ALLELE = 'ALLELE1',
                   OTHER_ALLELE = 'ALLELE0',
                   EAF = 'A1FREQ',
                   P = P_BOLT_LMM_INF) |>
            mutate(N = sample_size) |>
            left_join(HWE, by = "SNPID") |>
            left_join(info, by = "SNPID") |>
            select(all_of(columns))
  

dat |> write.table(output_file, quote = FALSE, 
                   col.names = TRUE, row.names = FALSE, sep = "\t")  
  
  
}


format_bolt(bolt_file = "data/BOTH_stats",
            pheno_file = "data/pheno.txt",
            info_file = "data/snp.info",
            hwe_file = "data/MCTFR.hwe",
            output_file = paste0("output/MCTFR.BOTH.association-results.", Sys.Date(), ".txt"))

format_bolt(bolt_file = "data/HIGH_stats",
            pheno_file = "data/first_order_pheno.txt",
            info_file = "data/snp.info",
            hwe_file = "data/MCTFR.hwe",
            output_file = paste0("output/MCTFR.BOTH.association-results.", Sys.Date(), ".txt"))

format_bolt(bolt_file = "data/LOW_stats",
            pheno_file = "data/second_order_pheno.txt",
            info_file = "data/snp.info",
            hwe_file = "data/MCTFR.hwe",
            output_file = paste0("output/MCTFR.BOTH.association-results.", Sys.Date(), ".txt"))

