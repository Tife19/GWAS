#!/usr/bin/env Rscript

# Load necessary library
#if (!Package("qqman", quietly = TRUE)) {
  #install.packages("qqman")  # Install qqman package if not already installed
#}
library(qqman)

# Run PLINK command (this should be executed in the terminal, not in R)
#system("./plink --bfile quantfamdata --assoc --linear --out quantfamdata_assoc")

# Read the results from the PLINK output
res3 <- read.table("quantfamdata_assoc.assoc.linear", header = TRUE)

# Display the first few rows of the results
print(head(res3))

# Calculate genomic inflation factor (lambda)
chi <- qchisq(1 - res3$P, 1)
lambda <- median(chi) / 0.456
cat("Lambda value for population structure assessment:", lambda, "\n")  # Print lambda value

# Create a new data frame with relevant columns
gwas_results <- data.frame(SNP = res3$SNP, CHR = res3$CHR, BP = res3$BP, P = res3$P)

# Display the first few rows of the new data frame
print(head(gwas_results))
#To check the number of snps on a chromosome
as.data.frame(table(gwasResults$CHR))

# Generate QQ plot and save as JPEG
jpeg("qq_quantfamdata.jpeg")
qq(gwas_results$P, main= "QQ plot for Association testing", xlim= c(0, 7), ylim= c(0, 30), pch = 18, col = "coral4", cex = 1.9, las = 1)

#Add lambda value to the plot
text(x = 1, y = 28, labels = paste("λ =", round(lambda, 3)), col = "black", cex = 1.2)

#Define significant threshold
genome_wide_threshold <- 5e-8
suggestive_threshold <-1e-5

#Identify significant and suggestive snps
significant_snps <- gwas_results[gwas_results$P < genome_wide_threshold, ]
suggestive_snps <- gwas_results[gwas_results$P < suggestive_threshold & gwas_results$P >= genome_wide_threshold, ]

# Save significant SNPs to a file
write.table(significant_snps, file = "significant_snps.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)
cat("Significant SNPs have been saved to 'significant_snps.txt'\n")

#Generate Manhattan plot and save as JPEG
jpeg("Manhanttan_quantfamdata.jpeg")
manhattan(gwas_results, main = "Manhattan plot for Association testing", snp="SNP", chr="CHR", bp="BP", p="P", ylim=c(0, 28),
cex = 1.5, col = c("blue4", "orange3"))


dev.off()

cat ("Plots have been generated and saved as JPEG files")
