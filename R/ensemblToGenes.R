#' Takes a dataframe and a column (or rownames) that has Ensembl IDs and returns dataframe with a column (or rownames) of gene names.
#' @param data A dataframe
#' @param column The name of the column that has the Ensembl IDs. If IDs are stored as rownames, use column = "rownames"
#' @param type Whether the type of Ensembl IDs are transcript or gene
#' @param keep.ensembl.ids Whether to keep Ensembl IDs as a separate column
#' @param make.genes.unique Choose whether you want only unique gene names. Duplicated genes will be separated by a dash. If IDs are stored as rownames,
#' make.genes.unique must be TRUE
#' @return A dataframe with column now as Ensembl IDs
#' @examples
#' ensemblToGenes(data = counts.table, column = "Genes", type = "gene")
#' @export

ensemblToGenes <- function(data, ...) {
  UseMethod("ensemblToGenes")
}

#' @export
ensemblToGenes.tbl_df <- function(data = data,
                                  column = column,
                                  type = type,
                                  keep.ensembl.ids = FALSE,
                                  make.genes.unique = TRUE,
                                  ...) {
  
  type <- stringr::str_to_lower(type)
  if (type == "gene"){
    
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
  
  if (isTRUE(keep.ensembl.ids)) {
    data <- dplyr::relocate(data, external_gene_name, .after = !!as.name(column))
    data <- dplyr::rename(data, "Ensembl_ID" := !!column)
    data <- dplyr::rename(data, "Gene_Name" = external_gene_name)
    
  } else {
    data <- dplyr::relocate(data, external_gene_name, .before = !!as.name(column))
    data <- dplyr::select(data,-!!as.name(column))
    data <- dplyr::rename(data, "Gene_Name" = external_gene_name)
  }
  if (isTRUE(make.genes.unique)) {
    data <- dplyr::mutate(data, Gene_Name = make.unique(Gene_Name, sep = "_"))
  }
  return(data)
}

#' @export
ensemblToGenes.data.frame <- function(data = data,
                                      column = column,
                                      type = type,
                                      keep.ensembl.ids = FALSE,
                                      make.genes.unique = TRUE,
                                      ...) {
  
  if (identical(column, "rownames") & make.genes.unique == FALSE){
    stop("For data frames with rownames as genes, make.genes.unique must
         be set to TRUE")
  }
  
  if (identical(column, "rownames")) {
    rownames <- rownames(data)
    data <- dplyr::mutate(data, "ensembl_ids" = rownames)
    rownames(data) <- NULL
    column <- "ensembl_ids"
    check.rows = TRUE
  }
  
  type <- stringr::str_to_lower(type)
  if (type == "gene"){
    
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
  
  if (isTRUE(keep.ensembl.ids)) {
    data <- dplyr::relocate(data, external_gene_name, .after = !!as.name(column))
    data <- dplyr::rename(data,  "Ensembl_ID" := !!column)
    data <- dplyr::rename(data, "Gene_Name" = external_gene_name)
    
  } else {
    data <- dplyr::relocate(data, external_gene_name, .before = !!as.name(column))
    data <- dplyr::select(data, -!!as.name(column))
    data <- dplyr::rename(data, "Gene_Name" = external_gene_name)
  }
  if (isTRUE(make.genes.unique)) {
    data <- dplyr::mutate(data, Gene_Name = make.unique(Gene_Name, sep = "_"))
  }
  
  if (isTRUE(check.rows)){
    gene.names <- magrittr::extract2(data, "Gene_Name")
    data <- dplyr::select(data, -c(Gene_Name))
    rownames(data) <- gene.names
  }
  
  return(data)
}

#' @export
ensemblToGenes.matrix <- function(data = data,
                                  column = column,
                                  type = type,
                                  keep.ensembl.ids = FALSE,
                                  make.genes.unique = TRUE,
                                  ...) {
  
  if (identical(column, "rownames") & make.genes.unique == FALSE){
    stop("For data frames/matrices with rownames as genes, make.genes.unique must
         be set to TRUE")
  }
  
  rownames <- rownames(data)
  data <- tibble::as_tibble(data, rownames = "ensembl_ids")
  column <- "ensembl_ids"
  
  
  type <- stringr::str_to_lower(type)
  if (type == "gene"){
    
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
  
  if (isTRUE(keep.ensembl.ids)) {
    data <- dplyr::relocate(data, external_gene_name, .after = !!as.name(column))
    data <- dplyr::rename(data,  "Ensembl_ID" := !!column)
    data <- dplyr::rename(data, "Gene_Name" = external_gene_name)
    
  } else {
    data <- dplyr::relocate(data, external_gene_name, .before = !!as.name(column))
    data <- dplyr::select(data, -!!as.name(column))
    data <- dplyr::rename(data, "Gene_Name" = external_gene_name)
  }
  if (isTRUE(make.genes.unique)) {
    data <- dplyr::mutate(data, Gene_Name = make.unique(Gene_Name, sep = "_"))
  }
  
  gene.names <- magrittr::extract2(data, "Gene_Name")
  data <- dplyr::select(data, -c(Gene_Name))
  data <- as.matrix(data)
  rownames(data) <- gene.names
  
  return(data)
}


