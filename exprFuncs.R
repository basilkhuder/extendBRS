cpmGeneExpr <- function(counts, gene.col, goi) { 
  counts <- tibble::column_to_rownames(counts, var = gene.col)
  counts <- edgeR::cpm(counts, log = TRUE)
  counts <- tibble::as_tibble(counts, rownames = gene.col)
  counts <- dplyr::filter(counts,!!as.name(gene.col) %in% c(goi))
  counts <- tidyr::pivot_longer(counts,cols = -all_of(gene.col),names_to = "Samples",
                                values_to = "CPM")
  return(counts)
  
} 

