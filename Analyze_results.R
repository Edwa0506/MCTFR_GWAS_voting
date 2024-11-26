setwd("/panfs/jay/groups/17/leej5/edwa0506/voting_gwas")

devtools::install_github('kaustubhad/fastman', build_vignettes = TRUE)
library(fastman)
library(vroom)

# Load the data
low <- vroom::vroom(paste0("output/MCTFR.LOW.association-results.", Sys.Date(), ".txt"))
high <- vroom::vroom(paste0("output/MCTFR.HIGH.association-results.", Sys.Date(), ".txt"))
both <- vroom::vroom(paste0("output/MCTFR.BOTH.association-results.", Sys.Date(), ".txt"))

# Save QQ plot for low
pdf("plots/qq_low.pdf", width = 7.5, height = 7.5)
fastqq(low$P)
dev.off()

# Save QQ plot for high
pdf("plots/qq_high.pdf", width = 7.5, height = 7.5)
fastqq(high$P)
dev.off()

# Save QQ plot for both
pdf("plots/qq_both.pdf", width = 7.5, height = 7.5)
fastqq(both$P)
dev.off()

# Save Manhattan plot for low
pdf("plots/man_low.pdf", width = 7.5, height = 7.5)
fastman(low, bp = "BP_b37")
dev.off()

# Save Manhattan plot for high
pdf("plots/man_high.pdf", width = 7.5, height = 7.5)
fastman(high, bp = "BP_b37")
dev.off()

# Save Manhattan plot for both
pdf("plots/man_both.pdf", width = 7.5, height = 7.5)
fastman(both, bp = "BP_b37")
dev.off()
