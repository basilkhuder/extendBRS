#' Creates a fictitious Seurat object based upon a set number of genes and barcodes
#' @param num.genes Number of genes to use
#' @param num.barcodes Number of barcodes to use
#' @param probs.zero Probabiltiy of getting a zero count amongst the genes
#' @return A Seurat object
#' @examples
#' fictitiousSeurat(num.genes = 1000, num.barcodes = 1000, probs.zero = NULL)
#' @export

fictitiousSeurat <- function(num.genes, num.barcodes, probs.zero = NULL) { 
  download.file("https://basilkhuder.s3.us-east-2.amazonaws.com/ensembl_id_and_gene_names.txt", 
                  destfile = "ensembl_id_and_gene_names.txt")
  ensembl.genes <- readr::read_tsv("ensembl_id_and_gene_names.txt")
  ensembl.genes <- dplyr::slice_n(ensembl.genes, n = num.genes)
  ensembl.genes <- dplyr::pull(ensembl.genes, external_gene_name)
  
   download.file("https://basilkhuder.s3.us-east-2.amazonaws.com/10XGenomics_PBMC_Barcodes.tsv.gz", 
                  destfile = "barcodes.tsv.gz")
  
  barcodes <- readr::read_tsv("barcodes.tsv.gz", col_names = "barcodes") 
  barcodes <- dplyr::slice_n(barcodes, n = num.barcodes)
  barcodes <- dplyr::pull(barcodes, barcodes)
  
   download.file("https://basilkhuder.s3.us-east-2.amazonaws.com/sc_counts_probability_vector.Robj", 
                  destfile = "sc_counts_probability_vector.Robj")
  
  load("sc_counts_probability_vector.Robj")
  gene.counts <- sample.int(4000, num.genes * num.barcodes, prob = probs)
  mat <- Matrix::Matrix(gene.counts, nrow = num.genes, ncol = num.barcodes)
  mat <- CreateSeuratObject(mat)
  return(mat) 
} 

