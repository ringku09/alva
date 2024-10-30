
#' Expand existing data by adding compound
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `expand_data()` is used to add one or more compounds data with the existing data.
#' Its process new compound data and update both **expression** and **metadata**. You can add
#' compound(s) only from [open TG-GATEs](https://toxico.nibiohn.go.jp/english/) database. The
#' available compound can be found by calling the function [`compounds()`].
#'
#'
#' @param compds A character vector contains compound name or its abbreviation to add.
#' @param expr_data A matrix of gene expression data of given/existing compounds.
#' @param attr_data Metadata of given/existing compounds, usually a data frame or tibble object.
#' @inheritParams process_compaunds
#'
#' @return A list.
#'   * A matrix of expanded gene expression data of compounds, probes in row and sample (CEL-ID) in column.
#'   * A data frame of expanded metadata about sample and phenotypic information of compounds.
#' @export
#'
#' @examples
#' \donttest{
#' xxx <- expand_data( compds, expr_data, attr_data)
#' }
expand_data <- function(compds,
                        expr_data,
                        attr_data,
                        multicore = FALSE,
                        store = FALSE,
                        method = "rma",
                        output_dir = missing_arg(),
                        error_call = caller_env()) {
  error_msg <- error_indata(expr_data, attr_data)
  if (!is.null(error_msg)) {
    cli_abort(error_msg)
  }
  output_dir <- destination(output_dir)
  new_comp <- comp_abbr(compds)
  ext_comp <- unique(attr_data$compound_abbr)
  ext_idx <- new_comp %in% ext_comp
  if (any(ext_idx)) {
    cli_warn(c(
      "Some of the compound already exist in the data",
      "i" = "{style_bold(col_red(abbr2name(new_comp[ext_idx])))} {?is/are} already in the data,
      we ignore {?it/them} to add again."
    ))
  }
  new_dl <- new_comp[!ext_idx]
  data_type <- get_dtype(attr_data)
  log2FC <- unique(attr_data$fc)
  if (length(log2FC) > 1) {
    cli_abort(c("Data must be unique type.",
                "x" = "You have used both `log2FC` and `normal` data.",
                "i" = "Please use only one type either `log2FC` or `normal` data.")
              )
  }
  new_data <- get_TGGATEs(compds = new_dl, data_type = data_type, fc = log2FC,
                          method = method, multicore = multicore,
                          store = store, output_dir = output_dir)
  expr_mat <- cbind(expr_data, new_data$expression)
  dict_df <- rbind(attr_data, new_data$metadata)
  if (store) {
    expr_mat %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    dict_df %>%
      readr::write_csv(file = paste(output_dir, "attribute.csv", sep = "/"))
  }
  tgx_data <- list(expr_mat, dict_df)
  names(tgx_data) <- c("expression", "attribute")
  return(tgx_data)
}


#' Shorten data by removing compounds
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `shirink_data()` is used to remove one or more compounds data from the existing data.
#' The corresponding samples of query compounds remove from both gene expression and metadata. The
#' compound(s) you want to remove must contain in both expression and attribute data.
#'
#' @inheritParams expand_data
#'
#' @return A list.
#'   * A matrix of shorten gene expression data of compounds, probes in row and sample (CEL-ID) in column.
#'   * A data frame of shorten metadata about sample and phnotypic information of compounds.
#'
#' @export
#'
#' @examples
#' \donttest{
#' compds <- c("APAP", "MP", "hfgjh")
#' xxx <- shrink_data(compds, expr_data, attr_data)
#' }
shrink_data <- function(compds,
                        expr_data,
                        attr_data,
                        store = FALSE,
                        output_dir = missing_arg()) {
  error_msg <- error_indata(expr_data, attr_data)
  if (!is.null(error_msg)) {
    cli_abort(error_msg)
  }
  output_dir <- destination(output_dir)
  rm_comp <- comp_abbr(compds)
  ret_idx <- attr_data$compound_abbr %in% rm_comp
  dict_df <- attr_data[!ret_idx, ]
  act_bar <- match(dict_df$barcode, colnames(expr_data))
  expr_mat <- expr_data[, act_bar]
  if (store) {
    expr_mat %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    dict_df %>%
      readr::write_csv(file = paste(output_dir, "attribute.csv", sep = "/"))
  }
  tgx_data <- list(expr_mat, dict_df)
  names(tgx_data) <- c("expression", "attribute")
  return(tgx_data)
}


#' Clean and refresh data
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `cnr_data()` is used to update your given/existing gene expression and meta data by
#' adding or removing compound(s) data based on query compound(s).
#'
#' The function `cnr_data()` is used to remove, or add some compound data with the existing data.
#' The compound(s) you want to add must available in open
#' [open TG-GATEs](https://toxico.nibiohn.go.jp/english/) database. A list of available compounds
#' can be found by calling the function [compounds()]. However, the compound(s) you want to remove must
#' contain in both expression and attribute data.
#'
#' @param ... Character vectors of compounds name or abbreviation in groups to analyse. Each argument
#'   can either be a character vector, a list that could be a character vector,  or a list of
#'   character vectors. Missing
#' @inheritParams expand_data
#'
#' @return A list.
#'   * A matrix of gene update expression measure of given compounds, probes in row and sample (CEL-ID) in column.
#'   * A data frame of update attribute about sample information of given compounds.
#' @export
#'
#' @examples
#' output_dir <- "D:/project1_TGx/data_processing5"
#'comps_gr <- list(GHS = c("BBZ", "APAP", "CMA", "MP", "NFZ"),
#'                 Non_GHS = c("EME", "GMC", "GBC", "HCB", "PEN", "INAH", "PH"))
#'xxx <- cnr_data(
#'  comps_gr,
#'  expr_data = expr_data,
#'  attr_data = attr_data,
#'  organism = "rat",
#'  data_type = "R1LS",
#'  log2FC = TRUE,
#'  multicore = FALSE,
#'  output_dir = "D:/project1_TGx/data_processing5",  # tempdir()
#'  store = FALSE
#')
cnr_data <- function(...,
                     expr_data,
                     attr_data,
                     probes = NULL,
                     multicore = FALSE,
                     store = FALSE,
                     method = "rma",
                     output_dir = missing_arg(),
                     error_call = caller_env()) {
  comps_group <- handl_group(...)
  comps_gr <- lapply(comps_group, comp_abbr, error_call = error_call)
  error_msg <- error_indata(expr_data, attr_data)
  if (!is.null(error_msg)) {
    cli_abort(error_msg)
  }
  compds <- as.character(unlist(comps_gr))
  compds <- comp_abbr(compds)
  mis_comps <- compds[!(compds %in% attr_data$compound_abbr)]
  attr_comp <- unique(attr_data$compound_abbr)
  ext_comps <- attr_comp[!(attr_comp %in% compds)]
  if (!is_empty(ext_comps)) {
    cli_alert_warning(c("The input data for the compound{?s}
                        {style_bold(col_red(backtick(abbr2name(ext_comps))))}
                        has been removed from the analysis because {?this/these}
                        input compounds are not listed in the query list."),wrap = TRUE)
    ext_data <- shrink_data(
      ext_comps,
      expr_data,
      attr_data,
      output_dir = output_dir,
      store = FALSE
      )
    expr_data <- ext_data$expression
    attr_data <- ext_data$attribute
  }
  if (!is_empty(mis_comps)) {
    exp_data <- expand_data(
      mis_comps,
      expr_data,
      attr_data,
      multicore = multicore,
      output_dir = output_dir,
      store = store)
    expr_data <- exp_data$expression
    attr_data <- exp_data$attribute
  }
  attr_new <- attr_data %>%
    dplyr::mutate(compound_abbr = factor(.data$compound_abbr, levels = compds)) %>%
    dplyr::arrange(.data$compound_abbr)
  tgx_class <- vector(length = nrow(attr_new))
  for(i in 1:length(comps_gr)) {
    tgx_class[attr_new$compound_abbr %in% comps_gr[[i]]] <- names(comps_gr)[i]
  }
  tgx_class <- factor(tgx_class, levels = names(comps_gr))
  attr_new <- attr_new %>%
    dplyr::mutate(group = tgx_class, .after = barcode)
  #unique(comp_dict$compound_abbr)
  act_bar <- match(attr_new$barcode, colnames(expr_data))
  expr_new <- expr_data[, act_bar]
  if (is.null(probes)) {
    expr_new <- expr_new
  } else{
    true_probes <- probes %in% rownames(expr_new)
    if (any(!true_probes)) {
      false_probes <- probes[!true_probes]
      n_false <- length(false_probes)
      cli_alert_warning(c("Given {.emph {n_false}} probe{?s} ",
                          "{style_bold(col_red(backtick(false_probes)))} ",
                          "{?is/are} not found in the gene expression data."),
                        wrap = TRUE)
    }
    expr_new <- expr_new[probes[true_probes],]
  }
  if (store) {
    expr_new %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    attr_new %>%
      readr::write_csv(file = paste(output_dir, "attribute.csv", sep = "/"))
  }
  return(list(expr_data = expr_new, attr_data = attr_new))
}
