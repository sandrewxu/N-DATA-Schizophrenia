library(Matrix)
source("DNM-Gene-List.R")   # Step 1: Create a disease de novo list from Gene4Denovo
source("dnm_sub-creator.R") # Step 2: Using a known gene list, create a list of de
                            # novo genes and direct neighbors (with names in STRING
                            # network) and find mutability
source("NDATA-Functions.R") # Step 3: Run NDATA

### Step 1
dnmGeneList <- get_dnmGeneList("SCZ")
dnmGeneList$Gene <- as.character(dnmGeneList$Gene)

# Organize the dnm gene list. After running once, determine the genes that do not
# match up with PPI and change here. These are crucial to analysis
dnmGeneList$Gene[which(dnmGeneList$Gene == "ACKR1")] <- "DARC"
dnmGeneList$Gene[which(dnmGeneList$Gene == "ADGRE1")] <- "EMR1"
dnmGeneList$Gene[which(dnmGeneList$Gene == "ADGRF5")] <- "GPR116"
dnmGeneList$Gene[which(dnmGeneList$Gene == "ADGRL2")] <- "LPHN2"
dnmGeneList$Gene[which(dnmGeneList$Gene == "ATP5F1A")] <- "ATP5A1"
dnmGeneList$Gene[which(dnmGeneList$Gene == "CFAP43")] <- "WDR96"
dnmGeneList$Gene[which(dnmGeneList$Gene == "DNAAF4")] <- "DYX1C1"
dnmGeneList$Gene[which(dnmGeneList$Gene == "ELP1")] <- "IKBKAP"
dnmGeneList$Gene[which(dnmGeneList$Gene == "FCSK")] <- "FUK"
dnmGeneList$Gene[which(dnmGeneList$Gene == "FPGT-TNNI3K;TNNI3K")] <- "TNNI3K"
dnmGeneList$Gene[which(dnmGeneList$Gene == "JMJD7-PLA2G4B;PLA2G4B")] <- "PLA2G4B"
dnmGeneList$Gene[which(dnmGeneList$Gene == "MRTFB")] <- "MKL2"
dnmGeneList$Gene[which(dnmGeneList$Gene == "MTCL1")] <- "SOGA2"
dnmGeneList$Gene[which(dnmGeneList$Gene == "NSD2")] <- "WHSC1"
dnmGeneList$Gene[which(dnmGeneList$Gene == "NSD3")] <- "WHSC1L1"
dnmGeneList$Gene[which(dnmGeneList$Gene == "OR2A1;OR2A42")] <- "OR2A42"
dnmGeneList$Gene[which(dnmGeneList$Gene == "PCNX2")] <- "PCNXL2"
dnmGeneList$Gene[which(dnmGeneList$Gene == "SEPTIN1")] <- "SEPT1"
dnmGeneList$Gene[which(dnmGeneList$Gene == "VIRMA")] <- "KIAA1429"

### Step 2
# Known Gene List (244 Genes): http://www.nature.com/articles/s41398-021-01294-x
# Known Gene List (167 Genes): https://www.medrxiv.org/content/10.1101/2020.09.18.20192815v1
# known_gene_list <- read.csv("known_genes_multiomics.txt", sep = "\t", header = F, col.names = "Gene")
# known_gene_list <- read.csv("~/Downloads/2021-22/Zhao Lab/NDATA/Schizophrenia/STRING/Known Genes/MedRXIV.txt", sep = "\t", header = F, col.names = "Gene")
known_gene_list$Count <- "0"

# Edge Cutoff 400, 95th confidence interval
ppi <- readRDS("string_ppi_400_weighted.RDS")

# Create a "base" gene list, which is a combination of the known and dnm gene lists
base_dnm_list <- unique(rbind(dnmGeneList, known_gene_list))
base_dnm_list <- base_dnm_list[order(base_dnm_list$Gene), ]

dnm_sub <- get_dnm_sub(ppi, base_dnm_list, dnmGeneList)

### Step 3
ppi_inter <- ppi[dnm_sub$Gene, dnm_sub$Gene]
dnm_sub$initValue <- ifelse(dnm_sub$Gene%in%known_gene_list$Gene, 1, -1)

#Network information, w0: no network info, w1: degrees in the network, w2: equal weight
d <- colSums(ppi_inter)
neighbors<-apply(ppi_inter,1,function(x) which(x>0))
w0 <- rep(0,length(d))
w1 <- sqrt(d)
w2 <- rep(0.5,length(d))

N <- 3402 #Probands
mu <- dnm_sub$Mutability
Y <- dnm_sub$Count

h_est <- c()
tau1_est <- c()
tau0_est <- c()
gamma_est <- c()
pseudollk <- c()

Error_fdr <- c()
Power_fdr <- c()
FDP_fdr <- c()


icm_res <- icm_randsearch(neighbors,w=w1,W=NULL,Y,mu,N,
                        theta_0_init = c(-4,0,0),theta_1_init = 3,threshold = 1e-4,
                        initLabel = dnm_sub$initValue)

h_est[i] <- icm_res$theta_0_est[1]
tau1_est[i] <- icm_res$theta_0_est[2]
tau0_est[i] <- icm_res$theta_0_est[3]
gamma_est[i] <- icm_res$theta_1_est

Gibbs <- Gibbs.post(icm_res$theta_0_est,icm_res$theta_1_est,neighbors,w=w1,W=NULL,
                    Y=Y,X=rep(1,length(Y)),mu,N,
                    mcmc_samples=500,burnin=200,thin=1,initLabel = icm_res$State)

pseudollk[i] <- mean(Gibbs$llk[201:500]/P)
q_0 <- Gibbs$q_0
State_post <- Get_Post_FDR(dnm_sub$Gene_match,q_0)$FDR
Error_fdr[i] <- sum(State_post[which(S[60,]==-1)]<0.05)/sum((S[60,]==-1))
FDP_fdr[i] <- sum(State_post[which(S[60,]==-1)]<0.05)/sum(State_post<0.05)
Power_fdr[i] <- sum(State_post[which(S[60,]==1)]<0.05)/sum((S[60,]==1))

result<-data.frame(pseudollk,gamma_est,
                   h_est,tau1_est,
                   Error_fdr,Power_fdr,FDP_fdr)
