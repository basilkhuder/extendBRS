#' Takes a dataframe and column that has Ensembl IDs and returns dataframe with column as Gene IDs.
#' @param data A dataframe
#' @param column The name of the column that has the Ensembl IDs. 
#' @param type Whether the type of Ensembl IDs are transcript or gene
#' @return A dataframe with column now as Ensembl IDs
#' @examples
#' ensemblToGenes(data = counts.table, column = "Genes", type = "gene")
#' @export

ensemblToGenes <- function(data = data, 
                           column = column,
                           type = type) {
  
  if(type == "gene"){ 
    
    download.file("https://basilkhuder.s3.us-east-2.amazonaws.com/ensembl_id_and_gene_names.txt", 
                  destfile = "ensembl_id_and_gene_names.txt")
    ensembl <- readr::read_tsv("ensembl_id_and_gene_names.txt")
    data <- dplyr::filter(data, !!as.name(column) %in% ensembl$ensembl_gene_id)
    
  } else if (type == "transcript"){ 
    
    download.file("https://basilkhuder.s3.us-east-2.amazonaws.com/ensembl_transcript_id_and_gene_names.txt", 
                  destfile = "ensembl_transcript_id_and_gene_names.txt")
    ensembl <- readr::read_tsv("ensembl_transcript_id_and_gene_names.txt")
    data <- dplyr::filter(data, !!as.name(column) %in% ensembl$ensembl_transcript_id)
    
  } else { 
    stop("Type must be either transcript or gene")
  }
  
  colnames(ensembl)[1] <- column
  data <- dplyr::inner_join(data, ensembl, by = column)
  data <- dplyr::select(data, c(external_gene_name, everything()), -column)
  data <- dplyr::rename(data, !!column := external_gene_name)
  return(data)
} 

