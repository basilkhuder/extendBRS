#Takes a dataframe and column that has Ensembl IDs and returns dataframe with column as Gene IDs. If an Ensembl ID does not have an associated Gene ID, that ID row will 
#be filtered out. 

ensemblToGenes <- function(data = data, 
                           column = column) {
  
  download.file("https://basilkhuder.s3.us-east-2.amazonaws.com/ensembl_id_and_gene_names.txt", 
                destfile = "ensembl_id_and_gene_names.txt")
  
  ensembl <- readr::read_tsv("ensembl_id_and_gene_names.txt")
  data <- dplyr::filter(data, !!as.name(column) %in% ensembl$ensembl_gene_id)
  colnames(ensembl)[1] <- column
  data <- dplyr::full_join(data, ensembl, by = column)
  data <- dplyr::select(data, c(external_gene_name, everything()), -column)
  data <- dplyr::rename(data, !!column := external_gene_name)
  return(data)
} 
