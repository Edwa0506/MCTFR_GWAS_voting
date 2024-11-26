#!/bin/bash -l
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=edwa0506@umn.edu 

cd /panfs/jay/groups/17/leej5/edwa0506/voting_gwas

#wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/old/BOLT-LMM_v2.3.tar.gz

#tar -xzvf BOLT-LMM_v2.3.tar.gz

#./BOLT-LMM_v2.3/bolt -h
# More recent versions don't work on the cluster. 
# Likely due to the below update made in v2.3.1
# Updated BLAS library to Intel MKL 2018 Update 1.


./BOLT-LMM_v2.3/bolt \
    --lmm \
    --bed=data/genotypes_rsid.bed \
    --bim=data/genotypes_rsid.bim \
    --fam=data/genotypes_rsid.fam \
    --phenoFile=data/first_order_pheno.txt \
    --phenoCol=pheno \
    --LDscoresFile=BOLT-LMM_v2.3/tables/LDSCORE.1000G_EUR.tab.gz \
    --modelSnps=data/ld_pruned_snps_for_pca.prune.in \
    --statsFile=data/HIGH_stats \
    --verboseStats