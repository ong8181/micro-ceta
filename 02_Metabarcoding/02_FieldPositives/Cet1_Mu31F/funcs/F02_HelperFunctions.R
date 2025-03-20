####
#### F02. Figure Helper functions
#### 2018.12.17 Ushio
#### 2022.02.19 Ushio
####

# taxa name summarize function
taxa_name_summarize <- function(ps_object,
                                taxa_rank,
                                top_taxa_n = 10,
                                taxa_rank_names = colnames(tax_table(ps_object))){
  tax_df <- as.data.frame(tax_table(ps_object))
  if(is.null(tax_df$rep_tax)) tax_df$rep_tax <- "Undetermined"

  # Search Others and Undetermined taxa
  rep_tax_cond1 <- tax_df[,taxa_rank] == "" & !is.na(tax_df[,taxa_rank])
  tax_col1 <- which(colnames(tax_df) == taxa_rank)
  tax_col2 <- which(colnames(tax_df) == "species")
  rep_tax_cond2 <- apply(tax_df[,(tax_col1+1):tax_col2] == "", 1, sum) == (tax_col2 - tax_col1)
  
  # Replace taxa names
  tax_df[!rep_tax_cond1, "rep_tax"] <- as.character(tax_df[!rep_tax_cond1, taxa_rank])
  tax_df[rep_tax_cond1 & !rep_tax_cond2, "rep_tax"] <- "Others"
  
  # Re-import phyloseq object with revised tax_table
  ps_object2 <- phyloseq(otu_table(ps_object), sample_data(ps_object), tax_table(as.matrix(tax_df)))
  
  # Repalce low abundance taxa name with Others
  taxa_abundance_rank <- aggregate(taxa_sums(ps_object2), by = list(tax_table(ps_object2)[,"rep_tax"]), sum)
  taxa_abundance_rank <- taxa_abundance_rank[order(taxa_abundance_rank$x, decreasing = T),]
  taxa_top <- taxa_abundance_rank[1:top_taxa_n,]
  
  low_tax <- is.na(match(tax_table(ps_object2)[,"rep_tax"], as.character(taxa_top[,1])))
  tax_table(ps_object2)[low_tax,"rep_tax"] <- "Others"
  
  return(ps_object2)
}


