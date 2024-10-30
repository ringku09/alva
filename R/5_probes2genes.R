#' Find gene information from probe set
#'
#' The function `probes2genes()` is used to get gene name, symbol, entrezid id etc from probes id.
#'
#' @param affy_ids Affymetric probe IDs
#' @param organism Sample organism used for experiment ("human" or "rat")
#'
#' @return
#' A table of Ensembl gene symbol and other information
#' @export
#'
#' @examples
#' \dontrun{
#' affyids=c("1368121_at", "1372510_at", "1384331_at", "1386958_at", "1390208_at", "1398791_at")
#' probe2gene(affy_ids = affyids, organism = "rat")
#' }
probes2genes <- function(affy_ids, organism = "rat") {
  if (identical(organism, "rat")) {
    gene_tab <- suppressMessages(AnnotationDbi::select(
      x = rat2302.db::rat2302.db,
      keys = affy_ids ,
      columns = c("SYMBOL","ENTREZID", "ENSEMBL","GENENAME")
    ))
  } else if(identical(organism, "human")) {
    gene_tab <- suppressMessages(AnnotationDbi::select(
      x = hgu133plus2.db::hgu133plus2.db,
      keys = affy_ids ,
      columns = c("SYMBOL","ENTREZID", "ENSEMBL","GENENAME")
    ))
  }
  gene_tab <- gene_tab[!duplicated(gene_tab$PROBEID), ]
  return(gene_tab)
}
