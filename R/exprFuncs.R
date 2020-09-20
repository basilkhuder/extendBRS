#' Calculate CPM values for a data.frame of counts and return specific genes only
#' @param counts Counts as tibble (gene or transcript ID in a column) 
#' @param gene.col The name of the column that has the gene names or transcript IDs. If genes are rownames, use gene.col = "rownames"
#' @param goi Genes interested in finding CPM values for 
#' @param log Whether you want logCPM or not. Default is TRUE
#' @return A tibble with CPM for genes of interest
#' @examples
#' cpmGeneExpr(counts, gene.col = "Genes", goi = c("Gene_1","Gene_2"))
#' @export

cpmGeneExpr <- function(counts, ...) { 
  UseMethod("cpmGeneExpr")
  } 

cpmGeneExpr.tbl_df <- function(counts, gene.col, goi, log = TRUE) { 
  counts <- tibble::column_to_rownames(counts, var = gene.col)
  if (isTRUE(log)){ 
    counts <- edgeR::cpm(counts, log = TRUE) 
  } else { 
        counts <- edgeR::cpm(counts, log = FALSE)
  } 
  counts <- tibble::as_tibble(counts, rownames = gene.col)
  counts <- dplyr::filter(counts, !!as.name(gene.col) %in% goi)
  counts <- tidyr::pivot_longer(counts,cols = -all_of(gene.col), names_to = "Samples",
                                values_to = "CPM")
  return(counts)
} 

cpmGeneExpr.data.frame <- function(counts, gene.col, goi, log = TRUE) { 
  if (!identical(gene.col, "rownames")) { 
    rownames(counts) <- magrittr::extract2(counts, gene.col) 
    counts <- dplyr::select(counts, -gene.col)
    } 
  
  if (isTRUE(log)){ 
    counts <- edgeR::cpm(counts, log = TRUE) 
  } else { 
        counts <- edgeR::cpm(counts, log = FALSE)
  } 
  
  if (identical(gene.col, "rownames")) { 
    return(counts[goi,])
    } 
  
  counts <- tibble::as_tibble(counts, rownames = gene.col)
  counts <- dplyr::filter(counts, !!as.name(gene.col) %in% goi)
  counts <- tidyr::pivot_longer(counts,cols = -all_of(gene.col), names_to = "Samples",
                                values_to = "CPM")
  return(counts)
} 

cpmGeneExpr.matrix <- function(counts, gene.col, goi, log = TRUE) { 
  
  if (isTRUE(log)){ 
    counts <- edgeR::cpm(counts, log = TRUE) 
  } else { 
        counts <- edgeR::cpm(counts, log = FALSE)
  } 
  return(counts[goi,])
} 

filterByExprAcross <- function(counts.list, 
                               groups.list, 
                               min.counts = 10, 
                               calc.norm.factors = FALSE, ...) { 
  UseMethod("filterByExprAcross")
}

filterByExprAcross.list <- function(counts.list, ...) {
  
  types <- purrr::map_chr(data, ~class(.x)[[1]])
  if (all(types == "data.frame")) {
    filterByExprAcross.data.frame(counts.list, groups.list, 
                              min.counts = min.counts, 
                              calc.norm.factors = calc.norm.factors)
  }
  
  if(all(types == "tbl_df")) { 
    filterByExprAcross.tbl.df(counts.list, groups.list, 
                              min.counts = min.counts, 
                              calc.norm.factors = calc.norm.factors)
    }
  
}

filterByExprAcross.tbl_df <- function(counts.list, 
                               groups.list, 
                               min.counts = 10, 
                               calc.norm.factors = FALSE){ 
  
  gene.col <- purrr::map(counts.list, ~ colnames(.x)[1])
  med.lib.size <- purrr::imap(counts.list, ~ dplyr::select(.x, -gene.col[[.y]]))
  med.lib.size <- purrr::map(med.lib.size, colSums)
  med.lib.size <- purrr::map(med.lib.size, median)
  min.counts.cpm <- purrr::imap(med.lib.size, ~ min.counts / .x * 1e6)
  cpm.counts <- purrr::imap(counts.list, ~ tibble::column_to_rownames(.x, var = gene.col[[.y]]))
  cpm.counts <- purrr::map(cpm.counts, ~ edgeR::cpm(.x, log = FALSE)) 
  cpm.counts <- purrr::imap(cpm.counts, ~ tibble::as_tibble(.x, rownames = gene.col[[.y]]))
  min.group <- purrr::map(groups.list, ~ min(table(.x)))
  
  filtered.genes <- purrr::imap(cpm.counts, ~ summarize.func(.x, gene.col[[.y]], min.counts.cpm[[.y]])) 
  filtered.genes <- purrr::reduce(filtered.genes, ~ dplyr::full_join(.x,.y,by = "Genes"))
  filtered.genes <- dplyr::filter(filtered.genes, !(sum.x < min.group[[1]] | sum.y < min.group[[2]]))
  filtered.genes <- dplyr::pull(filtered.genes, Genes)
  
  counts.list <- purrr::imap(counts.list, ~ dplyr::filter(.x, !!as.name(gene.col[[.y]]) %in% filtered.genes))
  counts.list <- purrr::imap(counts.list, ~ tibble::column_to_rownames(.x, var = gene.col[[.y]]))
  counts.list <- purrr::imap(counts.list, ~ edgeR::DGEList(.x, group = groups.list[[.y]]))
  
  if(isTRUE(calc.norm.factors)){
    counts.list <- dplyr::map(counts.list, edgerR::calcNormFactors)
  }
  
  return(counts.list)
}

summarizeFunc <- function(counts, 
                          gene.col, 
                          min.counts){ 
  filtered.genes <- dplyr::summarize(counts, dplyr::across(-{{gene.col}}, ~.x >= min.counts)) 
  filtered.genes <- dplyr::mutate(filtered.genes, 
                                  Genes = magrittr::extract2(counts, gene.col),
                                  sum = rowSums(filtered.genes))
  filtered.genes <- dplyr::select(filtered.genes, c(Genes, sum))
  return(filtered.genes)
}

