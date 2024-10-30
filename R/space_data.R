
get_space <- function(...,
                      expr_data,
                      attr_data,
                      probes = NULL,
                      dose = NULL,
                      time = NULL,
                      multicore = FALSE,
                      store = FALSE,
                      output_dir = missing_arg(),
                      error_call = caller_env()) {
  if (!is.null(dose) & !any(dose %in% c("Low", "Middle", "High"))) {
    cli_abort(c("The input dose level is incorrect.",
                "x" = "Dose level ({style_bold(col_red(backtick(dose)))}) not match with the data.",
                "i" = "Please choose dose level either `Low`, `Middle` or `High`.")
              , call = error_call)
  }
  if (!is.null(time) & !any(time %in% c("3 hr", "6 hr", "9 hr", "24 hr", "4 day", "8 day", "15 day", "29 day"))) {
    cli_abort(c("The input time level is incorrect.",
                "x" = "Time level ({style_bold(col_red(backtick(time)))}) not match with the data.",
                "i" = "Please choose time level either `3 hr`, `9 hr`, `9 hr` or `24 hr`.")
              , call = error_call)
  }
  error_msg <- error_indata(expr_data, attr_data)
  if (!is.null(error_msg)) {
    cli_abort(error_msg)
  }
  comps_gr <- list2(...)
  if (length(comps_gr) == 1 && is_bare_list(comps_gr[[1]])) {
    comps_gr <- comps_gr[[1]]
  }
  if (is_missing(output_dir)) {
    output_dir <- getwd()
  }
  comps_gr <- lapply(comps_gr, comp_abbr, error_call = error_call)
  compds <- unlist(comps_gr)
  ck_data <- cnr_data(
    comps_gr,
    expr_data = expr_data,
    attr_data = attr_data,
    probes = probes,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  expr_data <- ck_data$expr_data
  comp_dict <- ck_data$attr_data
#
#   if (is.null(probes)) {
#     expr_data <- expr_data
#   } else{
#     true_probes <- probes %in% rownames(expr_data)
#     if (any(!true_probes)) {
#       false_probes <- probes[!true_probes]
#       n_false <- length(false_probes)
#       cli_alert_warning(c("Given {.emph {n_false}} probe{?s} ",
#                           "{style_bold(col_red(backtick(false_probes)))} ",
#                           "{?is/are} not found in the gene expression data."),
#                         wrap = TRUE)
#     }
#     expr_data <- expr_data[probes[true_probes],]
#   }
  lev_str <- expr_str(comps_gr, attr_data = comp_dict)
  design_mat <- design_matrix(lev_str)
  if (is.null(names(comps_gr))) {
    names(comps_gr) <- paste("Group", LETTERS[1:length(comps_gr)], sep = "_")
  }

  if (is.null(time)) {
    expr_sp <- expr_data
  } else {
    comp_dict <- comp_dict %>%
      dplyr::filter(time_level == {{time}})
    expr_sp <- expr_data[, colnames(expr_data) %in% comp_dict$barcode]
  }
  if (is.null(dose)) {
    expr_sp <- expr_sp
  } else {
    comp_dict <- comp_dict %>%
      dplyr::filter(dose_level == {{dose}})
    expr_sp <- expr_sp[, colnames(expr_sp) %in% comp_dict$barcode]
  }

  #  For average data
  # if (identical(space, "dose")) {
  #   if (average) {
  #     expr_sp <- Y %*% design_mat$qGama
  #     #colnames(expr_sp) <- level_names(..., ..., ...)
  #   }
  #   dose_dict <- comp_dict %>%
  #     dplyr::filter(dose_level == "High")
  #   expr_sp <- Y[, colnames(Y) %in% dose_dict$barcode]
  #  dim(expr_sp)
  # }
  #
  # if (identical(space, "time")) {
  #   if (average) {
  #     expr_sp <- Y %*% design_mat$qDelta
  #     #colnames(expr_sp) <- level_names(..., ..., ...)
  #   }
  #   time_dict <- comp_dict %>%
  #     dplyr::filter(time_level == "24 hr")
  #   expr_sp <- Y[, colnames(Y) %in% time_dict$barcode]
  #   dim(expr_sp)
  # }
  tgx_class <- vector(length = nrow(comp_dict))
  for(i in 1:length(comps_gr)) {
    tgx_class[comp_dict$compound_abbr %in% comps_gr[[i]]] <- names(comps_gr)[i]
  }
  tgx_class <- factor(tgx_class, levels = names(comps_gr))
  comp_dict <- comp_dict %>%
    dplyr::mutate(group = tgx_class, .after = barcode)
  return(list(expr_data = expr_sp, attr_data = comp_dict))

}

# add class level in return



# probes <- res4$probe_id


avg_space <- function(...,
                 expr_data = expr_data,
                 attr_data =attr_data,
                 probes = NULL,
                 dose = NULL,
                 time = NULL) {
  space_data <- get_space(comps_gr,
                  expr_data = expr_data,
                  attr_data = attr_data,
                  probes = probes,
                  dose = dose,
                  time = time)
  space_nest <- space_data$attr_data %>%
    dplyr::group_by(compound_abbr, dose_level,time_level) %>%
    nest()
  space_barcd <- lapply(space_nest$data, function(x) x$barcode)
  barcd_expr <- lapply(space_barcd, function(x) rowMeans(as.matrix(space_data$expr_data[,x])))
  space_expr <- do.call(cbind, barcd_expr)

  doses <- substr(space_nest$dose_level, start=1, stop=1)
  timepoints <- gsub('([0-9]+) ([a-zA-Z])[a-zA-Z]*', '\\1\\2', space_nest$time_level)

  if (is.null(dose) & is.null(time)) {
    lab <- glue("{space_nest$compound_abbr} ({doses}-{timepoints})")
  } else if (is.null(dose)) {
    lab <- glue("{space_nest$compound_abbr} ({doses})")
  } else if (is.null(time)) {
    lab <- glue("{space_nest$compound_abbr} ({timepoints})")
  } else {
    lab <- space_nest$compound_abbr
  }
  space_attr <- space_nest %>%
    select(-data) %>%
    ungroup() %>%
    mutate(sample_id = lab)
  colnames(space_expr) <- lab
  return(list(expr_data = space_expr, attr_data = space_attr))
}
