###################################################################
###################################################################
###################################################################
###            Functions fgsea -> Created by Adrià           ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:    Working on it                                ######
###################################################################
###################################################################


if (!require("pacman", quietly = TRUE))
    install.packages("pacman")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fgsea")
library(fgsea)
pacman::p_load(dplyr,#### Analysis package
               ggplot2 #### Plot package
               ) 


GSEA = function(Deseq_results, GMT_data, p_value) {
  set.seed(1234)

  ##### Initial checks for fgsa
  if ( any( duplicated(names(Deseq_results)) )  ) {
    warning("Duplicates in gene names")
    Deseq_results = Deseq_results[!duplicated(names(Deseq_results))]
  }
  if  ( !all( order(Deseq_results, decreasing = TRUE) == 1:length(Deseq_results)) ){
    warning("Gene list not sorted")
    Deseq_results = sort(Deseq_results, decreasing = TRUE)
  }
  
  ### Run the fgsea
  fgRes <- fgsea::fgsea(pathways = GMT_data,
                        stats = Deseq_results,
                        minSize=1, ## minimum gene set size
                        maxSize=400, ## maximum gene set size
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(pval < !!p_value) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  #### Perform the collapsing to remove False positive Go codes enriched
  message("Collapsing Go codes -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = GMT_data,
                                      stats = Deseq_results)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated") ###Remove all 0 values
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  colos = setNames(c("#00A087FF", "#3C5488FF"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0("Top 10 (Total GO codes: Up=", sum(fgRes$Enrichment == "Up-regulated"),", Down=",    sum(fgRes$Enrichment == "Down-regulated"), ")")) 
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}
