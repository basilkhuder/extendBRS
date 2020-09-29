#' Takes a dataframe with a column of gene names, data frame with rownames as genes, or a vector of genes and returns dataframe with NCBI IDs (formerly Entrez) 
#' @param data A dataframe or a vector
#' @param column The name of the column that has the gene names. If genes are in rows, use column = "rownames", if vector use column = "vector"
#' @param drop.na Remove rows of genes that do not have an NCBI ID
#' @return Dataframe with NCBI IDs in a column
#' @examples
#' genesToEntrez(data, column = "Gene_Names", drop.na = TRUE) 
#' @export

genesToEntrez <- function(data, ...) {
  UseMethod("genesToEntrez")
}

genesToEntrez.tbl_df <- function(data, 
                          column,
                          drop.na = FALSE) { 
  
  genes_file <- tempfile(pattern = "ncbi_genes", fileext = ".txt")
  download.file(url = "https://basilkhuder.s3.us-east-2.amazonaws.com/ncbi_gene_id_and_gene_names.txt",
                destfile = genes_file)
  ncbi_ids <- readr::read_tsv(genes_file, col_names = TRUE, col_types = list(col_character(), col_character()))
  colnames(ncbi_ids[2]) <- column
  data <- dplyr::right_join(ncbi_ids, data)
  
  if (isTRUE(drop.na)){ 
    data <- tidyr::drop_na(data, NCBI_ID)
  }
  
  return(data)
}

genesToEntrez.data.frame <- function(data,
                                     column,
                                     drop.na = TRUE) { 
  
  if (identical(column, "rownames") & drop.na == FALSE){
    stop("For data frames with rownames as genes, drop.na must
         be set to TRUE")
    }
  
  genes_file <- tempfile(pattern = "ncbi_genes", fileext = ".txt")
  download.file(url = "https://basilkhuder.s3.us-east-2.amazonaws.com/ncbi_gene_id_and_gene_names.txt",
                destfile = genes_file)
  ncbi_ids <- readr::read_tsv(genes_file, col_names = TRUE, col_types = list(col_character(), col_character()))
  
  if (identical(column, "rownames")) {
    rownames <- rownames(data)
    data <- dplyr::mutate(data, "Gene_Name" = rownames)
    rownames(data) <- NULL
    check.rows = TRUE
    colnames(ncbi_ids[2]) <- "Gene_Name"
  } else { 
    colnames(ncbi_ids[2]) <- column
  }
  
  data <- dplyr::right_join(ncbi_ids, data, by = "Gene_Name")
  
  if (isTRUE(drop.na)) { 
    data <- tidyr::drop_na(data, NCBI_ID)
    }
  
  if (isTRUE(check.rows)) {
    data <- as.data.frame(data)
    gene.names <- magrittr::extract2(data, "Gene_Name")
    data <- dplyr::select(data, -c(Gene_Name))
    rownames(data) <- gene.names
    return(data)
  }
  
  data <- as.data.frame(data)
  return(data)
}

genesToEntrez.character <- function(data,
                                    column,
                                    drop.na = FALSE) { 
  
  genes_file <- tempfile(pattern = "ncbi_genes", fileext = ".txt")
  download.file(url = "https://basilkhuder.s3.us-east-2.amazonaws.com/ncbi_gene_id_and_gene_names.txt",
                destfile = genes_file)
  ncbi_ids <- readr::read_tsv(genes_file, col_names = TRUE, col_types = list(col_character(), col_character()))
  data <- dplyr::tibble(Gene_Name = data)
  data <- dplyr::right_join(ncbi_ids, data)
  return(data)
}
