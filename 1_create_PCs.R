setwd("/home/leej5/edwa0506/voting_gwas")
library(tidyverse)
# Create PCs
# Get maximal unrelated set
# Assign groups then perform PCA in plink


# subsetting to related takes a long time. 
# get perfectly genotypes variants to increase speed and accuracy 
info_scores <- vroom::vroom("data/snp.info")

# High info snps for Bolt mixed model
#info_scores |> filter(MAF > 0.01,
#                      Rsq > 0.9963) |>
#               select(rsid) |>
#               write.table("data/high_info_snps_for_bolt.txt",
#                           sep = " ", co)

info_scores |> 
  filter(Genotyped == "Genotyped") |>
  select(rsid) |>
  write.table("data/genotyped_snps.txt", sep = " ",
              col.names = FALSE, row.names = FALSE, quote =FALSE)

# Calc HWE
system2("./plink",
        args = c("--bfile", "data/genotypes_rsid",
                 "--hardy", 
                 "--out", "data/MCTFR"),
        stdout = TRUE,
        stderr = TRUE)

# Use non-imputed SNPs to finding related individuals faster
system2("./plink",
        args = c("--bfile", "data/genotypes_rsid",
                 "--extract", "data/genotyped_snps.txt",
                 "--make-bed", 
                 "--out", "data/genotypes_rsid_notimputed"),
        stdout = TRUE,
        stderr = TRUE)

# Find related individuals (used to exclude from calculating PC eigenvectors)
system2("./plink",
        args = c("--bfile", "data/genotypes_rsid_notimputed",
                 "--rel-cutoff", "0.025",
                 "--out", "data/related_individuals"),
        stdout = TRUE,
        stderr = TRUE)


fam_file <- read_table("data/genotypes_rsid.fam",
                       col_names =c("FID", "IID", "M", "P", "S", "pheno"))

unrelated <- read_table("data/related_individuals.rel.id", 
                      col_names = c("FID", "IID"))

fam_file$group <- ifelse(fam_file$IID %in% unrelated$IID, "g1", "g2") 

fam_file |> select(FID, IID, group) |> 
  write.table("data/related_cluster_file.txt",
              sep = " ",
              quote = FALSE,
              col.names = FALSE,
              row.names = FALSE) 

#FID, IID, cluster name 

# Prune SNPs to make PCs better and faster
system2("./plink",
        args = c("--bfile", "data/genotypes_rsid",
                 "--maf", "0.01",
                 "--indep-pairwise", "50", "10", "0.1",
                 "--out", "data/ld_pruned_snps_for_pca"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile",  "data/genotypes_rsid",
                 "--extract","data/ld_pruned_snps_for_pca.prune.in",
                 "--make-bed",
                 "--out", "data/genotypes_rsid_ld_pruned_pca"),
        stdout = TRUE,
        stderr = TRUE)

# Calculate PCs
system2("./plink",
        args = c("--bfile",  "data/genotypes_rsid_ld_pruned_pca",
                 "--within", "data/related_cluster_file.txt",
                 "--pca",
                 "--pca-cluster-names", "g1",
                 "--out", "data/pca"),
        stdout = TRUE,
        stderr = TRUE)
