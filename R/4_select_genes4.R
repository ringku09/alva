
#' Multistage hierarchical analysis of variance (MHANOVA) model
#'
#' The function `mhanova()` is used perform 4-stages hierarchical analysis of variance model.
#'
#' @param Y Gene expression data matrix, probes in row and sample (CEL-ID) in column.
#' @param design_mat Design matrices for MHANOVA model.
#' @param mse Logical parameter for output if MSE or p-values. default is `FALSE`.
#'
#' @return
#' A vector of p-values for different stage (group, compound, dose and time) or MSE.
#' @export
#'
#' @examples
#' \dontrun{
#' mhanova(expr_data, design_mat, mse = FALSE)
#' }
hlm <- function(Y, design_mat, p_adjust = TRUE, error_call = caller_env()) {
  if (!inherits(design_mat, "tggates")) {
    cli_abort(c("{.arg design_mat} must be a object of class `tggates`.",
                "x" = "The class {style_bold(col_cyan(backtick(class(design_mat))))} of \\
                {.var lev_str} is not supported.",
                "i" = "Please provide {.var design_mat} as a `tggates` object.")
              , call = error_call)
  }
  XX <- as(design_mat$R, 'dgCMatrix')
  SSE <- sse(Y, XX)
# SSE <- (Y %**% design_mat$R * Y) %**%  matrix(1, nrow = nrow(Y), ncol = 1) # 2nd
# SSE <- rowSums(Y %**% design_mat$R * Y)                                    # 3rd
# SSE <- diag(Y %**% design_mat$R %**% t(Y))                                 # not efficient
# SSE <- diag(tcrossprod((Y %**% design_mat$R), Y))                         # not efficient
  MSE <- SSE / (design_mat$N - design_mat$d...)
  SSA <- ssa(Y,design_mat$QA) # Y %**% design_mat$QA %.% Y
  MSA <- SSA / (design_mat$a-1)
  F_stat <- MSA/MSE
  p <- stats::pf(F_stat, design_mat$a - 1, design_mat$N - design_mat$d..., lower.tail = FALSE)
  # For t-test
  #t_stat <- sqrt(MSA/MSE)
  #p <- 2*pt(t_stat, df= design_mat$N - design_mat$d..., lower.tail = FALSE)
  if (p_adjust) {
    p <- p.adjust(p, method ="bonferroni", n = length(p))
  }
  names(p) <- rownames(Y)
  return(p)
}


#' Differentially expressed genes
#'
#' The function `de_genes()` is used to find differentially expressed (DE) genes.
#'
#' The function `de_genes()` is used to identify statistically significant genes
#' which are differentially expressed in expression levels between two or more compound groups.
#' A gene is said to be statistically significant if the *p-value* in all (or, some you interested) stages
#' (`group`, `compound`, `dose` and `time`) are significant.
#'
#' @param expr_data Gene expression data matrix
#' @param comp_attr Attribute data of the compounds
#' @param pval_cut P-value threshold
#' @param test_cond A character string to select *p-value* column(s) for significance test.
#'   If `NULL`(the default), all stages (`group`,`compound`, `dose` and `time`) will use for test.
#'   You can select stage(s) or column(s) by a character string contain stage name and separated
#'   by `&`, for instance, you want to test `group` and `dose` stage then character string will
#'   be `group&dose`, if you want `group`, `compound` and `time` then the character string will be
#'   `group&compound&time`, for single stage you can simply use the name of the column like `group`
#' @param p.cutof P-value cut off
#' @inheritParams process_compaunds
#'
#' @return
#' A `data.frame` or `tibble` about genes result.
#' @export
#'
#' @examples
#' \dontrun{
#' expr_data <- make_matrix(select(yy,-probes), dplyr::pull(yy, probes))
#' data_dict <- xx
#' compds <- list(p = c("BBZ", "APAP", "CMA", "MP", "NFZ"),
#'                N = c("EME", "GMC", "GBC", "HCB", "PEN", "INAH", "PH"))
#' res <- genes_MHANOVA(comps_gr = compds,
#'                      expr_data = expr_data,
#'                      attr_data = data_dict,
#'                      test_lev = "all",
#'                      p.cutof = 1e-8,
#'                      .parallel = FALSE)
#' }
de_genes <- function(...,
                     expr_data,
                     attr_data,
                     pval_cut = 1e-3,
                     gr_diff = FALSE,
                     multicore = FALSE,
                     store = FALSE,
                     output_dir = missing_arg(),
                     error_call = caller_env()) {
  comps_group <- handl_group(...)
  if (!length(comps_gr)>1) {
    cli_abort(c("The number of compound groups must be greater than one.",
                "i" = "Please use multiple compound groups for comparison."),
              wrap =TRUE, call = error_call)
  }
  error_msg <- error_indata(expr_data, attr_data)
  if (!is.null(error_msg)) {
    cli_abort(error_msg)
  }
  if (store) {
    if (is_missing(output_dir)) {
      output_dir <- getwd()
    }
  }
  comps_gr <- lapply(comps_group, comp_abbr, error_call = error_call)
  ck_data <- cnr_data(
    comps_gr,
    expr_data = expr_data,
    attr_data = attr_data,
    multicore = multicore,
    store = store,
    output_dir = missing_arg(),
    error_call = error_call
  )
  expr_data <- ck_data$expr_data
  lev_str <- expr_str(comps_gr, attr_data = ck_data$attr_data)
  design_mat <- design_matrix(lev_str)
  # Start parallel computing (if needed)
  # if (is.logical(multicore)) {
  #   if (multicore ) {
  #     multicore <- start_parallel(multicore)
  #     stop_cluster <- TRUE
  #   } else {
  #     multicore <- stop_cluster <- FALSE
  #   }
  # } else {
  #   stop_cluster <- if(inherits(multicore, "cluster")) FALSE else TRUE
  #   multicore <- start_parallel(multicore)
  # }
  # on.exit(if (multicore & stop_cluster)
  #   stop_parallel(attr(multicore, "cluster")))
  # # define operator to use depending on parallel being TRUE or FALSE
  # `%DO%` <- if (multicore) foreach::`%dopar%` else foreach::`%do%`
  # if (!multicore) {
  #   pval <- apply(expr_data, 1, mh_anova, design_mat = design_mat)
  #   pvalues <- t(pval)
  #   colnames(pvalues) <- c("group" ,"compound", "dose", "time")
  # } else {
  #   #i <- NULL
  #   ii <- attr(multicore, "cores")
  #   par_idx <- parallel::splitIndices(nrow(expr_data), ii)
  #   fun_call <- as.call(c(list(quote(foreach::foreach), i = seq_len(ii)), .combine ="cbind"))
  #   .fun <- eval(fun_call)
  #   pval <- .fun %DO% {apply(expr_data[par_idx[[i]], ], 1, mh_anova, design_mat = design_mat)}
  #   pvalues <- t(pval)
  #   colnames(pvalues) <- c("group" ,"compound", "dose", "time")
  # }
  pvalues <- hlm(Y = expr_data, design_mat = design_mat)
  #pmat <- select_pcol(pvalues, lab = test_cond)
  sig_probe <- pvalues <= pval_cut
  # Add code for average data
  #avg_gfc <- expr_data %**% design_mat$qAlfa
  sig_exp <- signif(pvalues, digits = 3)
  aov_df <- tibble::tibble(probe_id = names(sig_exp), p_values = sig_exp) %>%
    dplyr::mutate(sig_type = ifelse(sig_probe,"DE", "EE"))
  avg_gfc <- expr_data %**% design_mat$qAlfa
  colnames(avg_gfc) <- names(comps_gr)
  rownames(avg_gfc) <- rownames(expr_data)
  if (gr_diff) {
    aov_df <- col_diff(avg_gfc) %>%
      dplyr::left_join(aov_df, by = "probe_id")
  }
    aov_df <- tibble::as_tibble(avg_gfc, rownames = "probe_id") %>%
      dplyr::left_join(aov_df, by = "probe_id")
  return(aov_df)
}


#' Identify Toxigogenomics  Biomarker genes
#'
#' The function `tgx_genes()` is used to find final gene set using Leave One Out Compound (LOOC) method.
#' Genes misrepresented by a compound were removed using Leave One Out Compound (LOOC) method.
#' In LOOC process, each compound is removed from the experiment and find differentially expressed (DE)
#' genes for every LOO. The final gene set is therefore the intersection of the LOOC genes. Genes which
#' are co-regulated by a specific compound then filter out in LOOC process.
#'
#' @inheritParams de_genes
#'
#' @return
#' A data frame of gene identification result.
#' @export
#'
#' @examples
#' \dontrun{
#' expr_data <- make_matrix(select(yy,-probes), dplyr::pull(yy, probes))
#' data_dict <- xx
#' compds <- list(p = c("BBZ", "APAP", "CMA", "MP", "NFZ"),
#'                N = c("EME", "GMC", "GBC", "HCB", "PEN", "INAH", "PH"))
#' res <- LOO_MHANOVA(  comps_gr = compds,
#'                      expr_data = expr_data,
#'                      attr_data = data_dict,
#'                      test_lev = "all",
#'                      p.cutof = 1e-8,
#'                      .parallel = FALSE)
#' }
tgx_genes <- function(...,
                      expr_data,
                      attr_data,
                      pval_cut = 1e-3,
                      gr_diff = TRUE,
                      log_pvalue = TRUE,
                      multicore = FALSE,
                      store = FALSE,
                      output_dir = missing_arg(),
                      error_call = caller_env()) {
  comps_group <- handl_group(...)
  error_msg <- error_indata(expr_data, attr_data)
  if (!is.null(error_msg)) {
    cli_abort(error_msg)
  }
  if (store) {
    if (is_missing(output_dir)) {
      output_dir <- getwd()
    }
  }
  comps_gr <- lapply(comps_group, comp_abbr, error_call = error_call)
  check_data <- cnr_data(
    comps_gr,
    expr_data = expr_data,
    attr_data = attr_data,
    multicore = multicore,
    store = store,
    output_dir = missing_arg(),
    error_call = error_call
  )
  expr_data <- check_data$expr_data
  attr_data <- check_data$attr_data
  compds <- unlist(comps_gr)
  full_anova <- de_genes(
    comps_gr,
    expr_data = expr_data,
    attr_data = attr_data,
    pval_cut = pval_cut,
    gr_diff = gr_diff,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  chip <- unique(attr_data$arr_design)
  if(identical(chip, "Rat230_2")) {
    organism = "rat"
  }else if(identical(chip, "HG-U133_Plus_2")) {
    organism = "human"
  }
  aov_genes <- probes2genes(full_anova$probe_id, organism)
  # logFC <- col_diff(comps_gr,
  #                   expr_data[match(aov_genes$PROBEID, full_anova$probe_id), ],
  #                   design_mat$qAlfa)
  aov_tab <- full_anova[match(aov_genes$PROBEID, full_anova$probe_id), ]
  aov_df <- aov_tab %>%
    # dplyr::left_join(aov_tab, by = "probe_id") %>%
    dplyr::mutate(entrez_id = aov_genes$ENTREZID,
                  gene_name = aov_genes$SYMBOL,
                  gene_fnane = aov_genes$GENENAME,
                  .before = 1)
    # dplyr::relocate(sig_type, .after = probe_id)


  # arng_prob <- aov_tab %>%
  #   dplyr::select(c(.data$probe_id, .data$group, .data$compound, .data$dose, .data$time)) %>%
  #   tidyr::gather(label,  p_value,  c("group", "compound", "dose", "time")) %>%
  #   dplyr::group_by(probe_id) %>%
  #   dplyr::slice(which.min(p_value)) %>%
  #   arrange(p_value)
  # aov_tab2 <- aov_df[match(arng_prob$probe_id, aov_df$probe_id), ]
  aov_df <- aov_df %>%
    #na.omit() %>%
    distinct(gene_name, .keep_all= TRUE) %>%
    #distinct(probe_id, .keep_all= TRUE) %>%
    arrange(p_values)
  sig_probes <- aov_df$probe_id[aov_df$sig_type == "DE"]
  expr_loodata <- expr_data[sig_probes, ]
  loo_probes <- vector(mode = "list", length = length(compds))
  for (i in seq_len(length(compds))) {
    comp_grnew <- lapply(comps_gr, function(x) {
      if (compds[i] %in% x) x[x == compds[i]] else x })
    comp_loo <- attr_data$compound_abbr %in% unlist(comp_grnew)
    attr_new <- attr_data[comp_loo, ]
    expr_new <- expr_loodata[,match(attr_data$barcode[comp_loo], colnames(expr_loodata))]
    loo_anova <- de_genes(  # .parallel setup run each lapply , that take more time,
      comp_grnew,    #   need to setup once
      expr_data = expr_new,
      attr_data = attr_new,
      pval_cut = pval_cut,
      gr_diff = FALSE,
      multicore = FALSE,
      store = FALSE,
      output_dir = output_dir,
      error_call = error_call
    )
    loo_probes[[i]] <- loo_anova$probe_id[loo_anova$sig_type == "DE"]
  }
  comon_probes <- Reduce(intersect, loo_probes)
  sig_tab <- aov_df %>%
    dplyr::mutate(sig_type = ifelse(!(aov_df$probe_id %in% comon_probes) & (aov_df$sig_type == "DE"),
                                    "CE", sig_type)) %>%
    dplyr::relocate(probe_id, .before = entrez_id)
  # sig_tab <- loo_df
  # sig_tab <- proc_gtab(loo_df, full_p.cutof) %>%
  #   arrange(p_value)
  if (log_pvalue) {
    sig_tab <- sig_tab %>%
       dplyr::mutate(log10p = -log10(p_values),.after = p_values)
    # sig_tab <- loo_df %>%
    #   dplyr::mutate(dplyr::across(c(group, compound, dose, time), function(x) -log10(x)))
  }
  n_probe <- nrow(expr_data)
  aov_n <- sum(sig_tab$sig_type == "DE" | sig_tab$sig_type == "CE")
  loo_n <- sum(sig_tab$sig_type == "DE")
  fil1 <- n_probe-aov_n
  fil2 <- aov_n-loo_n
  # if (is.null(top)) {
    #cli::cli_alert_info("Identify {loo_n} genes at p-value threshold {pval_cut}")
    cli_ul(c(
      bg_blue("Total = {n_probe} probes"),
      bg_cyan("Step-1 filter = {fil1} probes"),
      bg_cyan("Step-2 filter = {fil2} probes"),
      bg_magenta("Identified = {loo_n} genes")
    ))
    #    cli::cli_alert_info("The LOO-MHANOVA identify {loo_n} genes at p-value threshold {loo_p.cutof}")
  #} else {
  #  sig_tab <- sig_tab %>%
  #    dplyr::mutate(sig_type = ifelse(
  #      (.data$probe_id %in% .data$probe_id[1:top]) & (.data$sig_type == "S"),
   #     "HS", sig_type))
   #     cli::cli_alert_info("The MHANOVA identify {aov_n} genes at p-value threshold {pval_cut}")
   # cli::cli_alert_info("Identify {loo_n} genes at p-value threshold {p.cutof} and top {top} genes
  #                      is select at p-value threshold {sig_tab$group[top]}")
  #}
  return(sig_tab)
}
