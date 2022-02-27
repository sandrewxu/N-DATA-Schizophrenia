get_dnm_sub <- function(ppi, base_dnm_list, dnmGeneList) {

  # 1. Match genes to row numbers in STRINGdb, find direct neighbors, and filter
  # for a 95th percentile confidence interval. Then, find the intersection with ppi.
  Candidate_id <- which(rownames(ppi) %in% base_dnm_list$Gene)
  dnm_sub <- unique(unlist(lapply(Candidate_id, function(x)
    row.names(ppi)[which(ppi[x,]>950)])))

  dnm_sub <- setdiff(dnm_sub, dnmGeneList$Gene)
  dnm_sub <- data.frame("Gene" = dnm_sub, "Count" = 0)

  ppi_genes <- unique(intersect(dnm_sub$Gene, row.names(ppi)))
  dnm_sub <- dnm_sub[which(dnm_sub$Gene%in%ppi_genes),]

  dnm_sub <- dnm_sub[order(dnm_sub$Gene),]
  dnm_sub <- dnm_sub[-1,]

  # 2. Match mutabilities
  mutability_table <- read.csv("EncoreDNM_mutability_table.txt", sep = "\t")
  mutability_table <- data.frame("Gene" = mutability_table$gene, "Mutability" = mutability_table$lof + mutability_table$Dmis)

  mutability_table$Gene[which(mutability_table$Gene == "ACKR1")] <- "DARC"
  mutability_table$Gene[which(mutability_table$Gene == "ADGRE1")] <- "EMR1"
  mutability_table$Gene[which(mutability_table$Gene == "ADGRF5")] <- "GPR116"
  mutability_table$Gene[which(mutability_table$Gene == "ADGRL2")] <- "LPHN2"
  mutability_table$Gene[which(mutability_table$Gene == "CFAP43")] <- "WDR96"
  mutability_table$Gene[which(mutability_table$Gene == "ELP1")] <- "IKBKAP"
  mutability_table$Gene[which(mutability_table$Gene == "MTCL1")] <- "SOGA2"
  mutability_table$Gene[which(mutability_table$Gene == "NSD2")] <- "WHSC1"
  mutability_table$Gene[which(mutability_table$Gene == "NSD3")] <- "WHSC1L1"
  mutability_table$Gene[which(mutability_table$Gene == "PCNX2")] <- "PCNXL2"
  mutability_table$Gene[which(mutability_table$Gene == "TENM1")] <- "TEN1"
  mutability_table$Gene[which(mutability_table$Gene == "VIRMA")] <- "KIAA1429"

  mutability_table <- mutability_table[order(mutability_table$Gene),]

  ppi_genes <- intersect(dnm_sub$Gene, mutability_table$Gene)
  dnm_sub <- dnm_sub[which(dnm_sub$Gene%in%ppi_genes),]

  dnm_sub$Mutability = mutability_table$Mutability[which(mutability_table$Gene%in%dnm_sub$Gene)]

  dnm_sub
}
