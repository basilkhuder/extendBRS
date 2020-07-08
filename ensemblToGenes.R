ensemblToGenes <- function(dataframe = data, 
                           column = column) {
  
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  ensembl <- biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                            filters = 'ensembl_gene_id',
                            values = data, 
                            mart = ensembl)
                            
  counts <- dplyr::filter(counts, Chr %in% ensembl$ensembl_gene_id)
  colnames(ensembl)[1] <- column
  counts <- full_join(counts, ensembl, by = column)
  counts <- dplyr::select(counts, c(external_gene_name, everything()), -column)
  counts <- dplyr::rename(counts, !!column := external_gene_name)
  return(counts)
} 
