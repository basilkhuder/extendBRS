#Takes a dataframe and column that has Ensembl IDs and returns dataframe with column as Gene IDs. If an Ensembl ID does not have an associated Gene ID, that ID row will 
#be filtered out. 

ensemblToGenes <- function(data = data, 
                           column = column) {
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  ensembl <- biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                            filters = 'ensembl_gene_id',
                            values = data[,column], 
                            mart = ensembl)
  data <- dplyr::filter(data, !!as.name(column) %in% ensembl$ensembl_gene_id)
  colnames(ensembl)[1] <- column
  data <- full_join(data, ensembl, by = column)
  data <- dplyr::select(data, c(external_gene_name, everything()), -column)
  data <- dplyr::rename(data, !!column := external_gene_name)
  return(data)
} 
