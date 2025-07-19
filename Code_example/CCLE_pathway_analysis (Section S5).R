library(clusterProfiler)
library(org.Hs.eg.db)

top10_factor1 <- c("UGT1A6", "FABP5", "PI3", "CXCL1", "IFITM1", 
                   "KRT17", "CYP4F11", "CXCL5", "SERPINB2", "DHRS2")

top10_factor2 <- c("CLU", "PI3", "RAC2", "KRT7", "ELF3", 
                   "CDKN2A", "EFEMP1", "TFPI2", "RPS4Y1", "SRGN")

top10_combined <- c("PI3", "FABP5", "CLU", "UGT1A6", "IFITM1", 
                    "TNFRSF6B", "S100A9", "CXCL1", "EFEMP1", "CXCL5")


# Top 13 genes
top13_factor1 <- c("UGT1A6", "FABP5", "PI3", "CXCL1", "IFITM1", "KRT17", "CYP4F11", 
                   "CXCL5", "SERPINB2", "DHRS2", "CD74", "NNMT", "PRSS3")

top13_factor2 <- c("CLU", "PI3", "RAC2", "KRT7", "ELF3", "CDKN2A", "EFEMP1", 
                   "TFPI2", "RPS4Y1", "SRGN", "TNFRSF6B", "XAGE1B", "ALDH3A1")

# Top 15 genes
top15_factor1 <- c("UGT1A6", "FABP5", "PI3", "CXCL1", "IFITM1", "KRT17", "CYP4F11", 
                   "CXCL5", "SERPINB2", "DHRS2", "CD74", "NNMT", "PRSS3", "TNFRSF6B", "S100A9")

top15_factor2 <- c("CLU", "PI3", "RAC2", "KRT7", "ELF3", "CDKN2A", "EFEMP1", 
                   "TFPI2", "RPS4Y1", "SRGN", "TNFRSF6B", "XAGE1B", "ALDH3A1", "S100A9", "TACSTD2")

# Top 20 genes
top20_factor1 <- c("UGT1A6", "FABP5", "PI3", "CXCL1", "IFITM1", "KRT17", "CYP4F11", 
                   "CXCL5", "SERPINB2", "DHRS2", "CD74", "NNMT", "PRSS3", "TNFRSF6B", "S100A9", 
                   "MAGEA4", "FN1", "ANKRD1", "BST2", "EFEMP1")

top20_factor2 <- c("CLU", "PI3", "RAC2", "KRT7", "ELF3", "CDKN2A", "EFEMP1", 
                   "TFPI2", "RPS4Y1", "SRGN", "TNFRSF6B", "XAGE1B", "ALDH3A1", "S100A9", "TACSTD2",
                   "KYNU", "FABP5", "IFITM1", "CPLX2", "FGB")


convert_to_entrez <- function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}

genes_factor1_entrez <- convert_to_entrez(top10_factor1)
genes_factor2_entrez <- convert_to_entrez(top10_factor2)
genes_combined_entrez <- convert_to_entrez(top10_combined)

kegg_factor1 <- enrichKEGG(gene = genes_factor1_entrez, organism = 'hsa', pAdjustMethod = "BH")
kegg_factor2 <- enrichKEGG(gene = genes_factor2_entrez, organism = 'hsa', pAdjustMethod = "BH")
kegg_combined <- enrichKEGG(gene = genes_combined_entrez, organism = 'hsa', pAdjustMethod = "BH")

print(summary(kegg_factor1))
print(summary(kegg_factor2))
print(summary(kegg_combined))
