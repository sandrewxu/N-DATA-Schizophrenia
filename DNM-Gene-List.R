### Pipeline Step 1: Take an De Novo Gene database and sort it to fit our needs
### Import Gene Database (Made for Gene4Denovo)
### Information: http://www.genemed.tech/gene4denovo/source
### Source: http://www.genemed.tech/gene4denovo/download
### 'De novo mutations and annotations in coding region', Version 1.1

get_dnmGeneList <- function (disease) {
  # Sort by disease
  geneDatabase <- read.csv("De_novo_mutations_and_annotations_in_coding_region_1.2.txt", sep = "\t", header = TRUE)
  geneDatabase <- geneDatabase[geneDatabase$Phenotype == disease,]

  # Add a column MutationClass, which is the primary filter
  func <- geneDatabase$Func.refGene
  exonicfunc <- geneDatabase$ExonicFunc.refGene
  MetaSVM <- substr(geneDatabase$MetaSVM, 1, 5)

  geneDatabase$MutationClass <- "NA"
  for (i in 1:length(geneDatabase$MutationClass)) {
    if (func[i] == "splicing") {geneDatabase$MutationClass[i] = "LoF"}
    if (func[i] == "exonic" || func[i] == "exonic;splicing") {
      if (exonicfunc[i] == "frameshift deletion" ||
          exonicfunc[i] == "frameshift insertion" || exonicfunc[i] == "stopgain" ||
          exonicfunc[i] == "stoploss") {geneDatabase$MutationClass[i] = "LoF"}
      else if (exonicfunc[i] == "synonymous SNV") {geneDatabase$MutationClass[i] = "Syn"}
      else if (exonicfunc[i] == "nonsynonymous SNV") {
        if (MetaSVM[i] > 0) geneDatabase$MutationClass[i] = "Dmis"
        else geneDatabase$MutationClass[i] = "Tmis"
      }
      else geneDatabase$MutationClass[i] = "NA"
    }
  }

  # Keep only the Dmis and LoF Mutations for Analysis
  # No need for mutability currently (will need to add before running NDATA)
  geneDatabase <- subset(geneDatabase, MutationClass == "Dmis" | MutationClass == "LoF")
  dnmGenes <- data.frame(table(geneDatabase$Gene.refGene))
  colnames(dnmGenes) <- c("Gene", "Count")

  dnmGenes
}
