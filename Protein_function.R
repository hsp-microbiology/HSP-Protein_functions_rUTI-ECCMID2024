###################################################################
###################################################################
###################################################################
###     Functions Main file exam -> Created by Adrià         ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:    Working on it                                ######
###################################################################
###################################################################

rm(list=ls())

############
#### Protein diferential analysis
############
Deseq_fecal <- readRDS() ### Read the samples, the matrix has to be Proteins on row and Samples on columns
Sample_data <- read.csv() ### Put the metadata of  samples, rows SAmples, columns variables of interest

##### Chech names
table(colnames(Deseq_fecal) == rownames(Sample_data))  ## All true

###Create de Deseq object

dds <- DESeqDataSetFromMatrix(countData = Deseq_fecal,
                              colData = Sample_data,
                              design = ~ condition)

# perform a minimal pre-filtering to remove rows that have only 0 or 1 read. 
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds , test = "Wald")
res <- results(dds)
resSig <- subset(res, padj < 0.05) ### Select all proteins P adjusted significant 

### Write the table for fecal sample::
write.csv(res,file = "Fecal_DESEQ.csv") ###Save the run


############
#### Protein enrichment analysis
############
# load in GO function info. from a csv file and get all de Go codes from samples, we need the proteins in dataset (minimum the selected by deseq2) and the GO codes macthed with that proteins.
pathData <- readRDS("Go_terms.rds")  ##Go terms, functions  

####Create a list with the proteins in each Go code
GMT_data <- list()

for (i in Go_codes) {
  EC1_go <- str_split_fixed(i,"\\[",2)[,2]
  EC1_go <- str_split_fixed(EC1_go,"\\]",2)[,1]
  GMT_data[[length(GMT_data) + 1]] <- pathData$Uniprot_ID[grep(paste(EC1_go),pathData$Gene.Ontology..GO.)]
}

names(GMT_data) <- Go_codes
### Get the LFC and proteins names from Deseq and create a numeric vector
Deseq_results <- resSig$log2FoldChange
names(Deseq_results) <- resSig$X

#####Chargue the function
source("./Functions/Enrichment.R)

####Run fgsea
res1 = GSEA(Deseq_results, GMT_data, p_value = 0.05)

#### Plot the results (only on interactive view)
res1[[2]] ###Save with ggsave()

#### Save GSEA results
results_GSEA <- as.data.frame(res1[[1]])
openxlsx::write.xlsx(results_GSEA,file = "Fecal_enrichment2.xlsx",rowNames = T)



