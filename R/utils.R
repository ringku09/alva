

#' Function for matrix of one
#'
#' The function `mat_one()` is used to generate a matrix of one.
#'
#' @param n Number of rows
#' @param m Number of columns
#'
#' @return A matrix of ones
#' @export
#'
#' @examples
#' mat_one(3, 4)
mat_one <- function(n, m) {
  matrix(1, nrow = n, ncol = m)
}

#' Direct sum of two vectors or matrices
#'
#' The function ``
#' @param A A numeric matrix or vector
#' @param B A numeric matrix or vector
#'
#' @return
#' a matrix of direct sum
#' @export
#'
#' @examples
#' A <- matrix(1, nrow = 3, ncol = 3)
#' B <- matrix(2, nrow = 3, ncol = 3)
#' direct_sum(A, B)
direct_sum <- function(A, B) {
  if (is.matrix(A) & is.matrix(B)) {
    A <- A
    B <- B
  } else if (is.vector(A) & is.vector(B)) {
    A <- as.matrix(A)
    B <- as.matrix(B)
  } else stop( "Argument A and B must be a matrix or vector" )
  up_mat <- cbind(A, matrix(0, nrow = nrow(A), ncol = ncol(B)))
  down_mat <- cbind(matrix(0, nrow = nrow(B), ncol = ncol(A)), B)
  dir_sum <- rbind(up_mat, down_mat)
  return(dir_sum)
}

# Set p-value sign based on coordinate
plot_pval <- function(label, p_value) {
  if (identical(label, "group")) {
    p_val <- p_value
  } else if (identical(label, "compound")) {
    p_val <- p_value
  } else if (identical(label, "dose")) {
    p_val <- -p_value
  } else if (identical(label, "time")) {
    p_val <- -p_value
  }
  return(p_val)
}

# Set log FC sign based on coordinate
plot_lfc <- function(label, logFC) {
  if (identical(label, "group")) {
    lfc <- logFC
  } else if (identical(label, "compound")) {
    lfc <- -logFC
  } else if (identical(label, "dose")) {
    lfc <- -logFC
  } else if (identical(label, "time")) {
    lfc <- logFC
  }
  return(lfc)
}

# Find middle potion of bar plot for group to display text
gr_barpos <- function(gr_data, var) {
  zz <- gr_data %>%
    select({{var}})
  z <- as.vector(unlist(zz))
  if(all(z >= 0)) {
    txt_lev <- rev(cumsum(lag(rev(z), default = 0)) + rev(z)/2)
  } else if (any(z <  0) & any(z >= 0)) {
    txt_lev <- vector("numeric", length(z))
    z1 <- rev(cumsum(lag(rev(z[z >= 0]), default = 0)) + rev(z[z >= 0])/2)
    z2 <- cumsum(lag(rev(z[z < 0]), default = 0)) + rev(z[z < 0])/2
    txt_lev[z >= 0] <- z1
    txt_lev[z < 0] <- z2
  } else if (all(z < 0)) {
    txt_lev <- vector("numeric", length(z))
    z1 <- rev(cumsum(lag(rev(z[z >= 0]), default = 0)) + rev(z[z >= 0])/2)
    z2 <- cumsum(lag(rev(z[z < 0]), default = 0)) + rev(z[z < 0])/2
    txt_lev[z >= 0] <- z1
    txt_lev[z < 0] <- z2
    txt_lev <- rev(txt_lev)
  }
  return(txt_lev)
}

# Dataframe to matrix
df2matrix <- function(df, row_names = NULL) {
  df_mat <-  as.matrix(df)
  if (!is.null(rownames))
    rownames(df_mat) = row_names
  return(df_mat)
}

# Average fold change
avg_fc <- function(x, y, FC = TRUE) {
  avg_gr <- x %*% y$qAlfa
  if (FC) {
    fc <- avg_gr[1] - avg_gr[2]
  } else fc <- avg_gr[1]
  return(fc)
}

# destination path creator
destination <- function(output_dir) {
  if (is_missing(output_dir)) {
    output_dir <- tempdir()
  }
  if(!file.exists(output_dir)) {
    dir.create(output_dir)
  }
  return(output_dir)
}

#' @title Compute pairwise difference between matrix columns
#' @param x A data matrix of size n times p. Where rows are observations and
#' columns are features.
#' @return A matrix of size n times (p choose 2), where each column is the
#' difference between two of the original columns.
#' @importFrom magrittr set_colnames
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr imap
#' @export
#' @examples
#' n = 10
#' p = 9
#' x = matrix(rep(1:p, n), nrow = n, ncol = p, byrow = TRUE)
#' rownames(x) <- paste0("r", 1:n)
#' y <- c("a", "b", "c")
#' names(y) <- paste0("X", 1:length(y))
#' dm <- cbind(c(1, 1, 1, 0, 0, 0, 0, 0, 0), c(0, 0, 0, 1, 1, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 0, 1, 1, 1))
#' col_diff(y, x, dm)
col_diff <- function(avg_data) {
  p <- ncol(avg_data)
  xin_list <- purrr::map(.x = seq_len(p - 1),
                         .f = ~ avg_data[, .x] - avg_data[, -c(seq_len(.x)), drop = FALSE])
  names(xin_list) <- colnames(avg_data)[-length(colnames(avg_data))]
  xin_list_name <- purrr::imap(
    .x = xin_list,
    .f = function(x, y) {
      mat_name = paste0(y, "-", colnames(x))
      colnames(x) = mat_name
      return(x)
    }
  )

  pair_diff <- do.call(cbind, xin_list_name)
  pair_mat <- tibble::as_tibble(pair_diff, rownames = "probe_id")
  return(pair_mat)
}


# with common legend
com_legend <- function(gg_plot) {
  tmp <- ggplot_gtable(ggplot_build(gg_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Capitalize first letter
block_fst <- function(x) {
  x <- tolower(x)
substr(x, 1, 1) <- toupper(substr(x, 1, 1))
x
}

# test expression and attribute data
error_indata <- function(expr_data, attr_data) {
  if (!inherits(expr_data, c("matrix", "array"))) {
    msg <- c("The expression data must be a matrix.",
                "x" = "The class {style_bold(col_cyan(backtick(class(expr_data))))} of \\
                expression data {?is/are} not supported.",
                "i" = "Please make expression data as matrix \\
             (`probes` in rows and `barcode` in columns).")

  }
  if (!inherits_any(attr_data, c("tbl_df", "tbl", "data.frame"))) {
    msg <- c("The attribute data must be a data frame.",
                "x" = "The class {style_bold(col_cyan(backtick(class(attr_data))))} of \\
                attribute data {?is/are} not supported.",
                "i" = "Please make attribute data as data frame.")

  }
  if (inherits(expr_data, c("matrix", "array")) &
      inherits_any(attr_data, c("tbl_df", "tbl", "data.frame"))) {
    match_bar <- colnames(expr_data) %in% attr_data$barcode
    if (any(!match_bar)) {
      msg <- c("All columns/barcodes in expression data must given in attribute data .",
               "x" = "{style_bold(col_cyan(backtick(colnames(expr_data)[match_bar])))} column{?s} \\
                {?is/are} not described in the attribute data.",
               "i" = "Please remove unspecified columns from the expression data. After that you \\
               can update your data by calling the function `cnr_data()`")
    } else msg <- NULL
  }
  return(msg)
}

is.compound <- function(comp_name, error_call = caller_env()) {
  if (is_empty(comp_name)) {
    cli_abort(c("{.var comp_name} must be non-empty.",
                "i" = "You have supplied an empty vector, please provide `name` or, `abbreviation` instred.")
              , call = error_call)
  }
  comp_tg <- compounds()
  comp_is <- comp_name %in% comp_tg$COMPOUND_NAME
  if (all(comp_is)) {
    idx <- match(comp_name, comp_tg$COMPOUND_NAME)
    abbr <- comp_tg$COMPOUND_ABBREVIATION[idx]
  }
  abbr_is <- comp_name %in% comp_tg$COMPOUND_ABBREVIATION
  if (all(abbr_is)) {
    abbr <- comp_name
  }
  if (all(!comp_is) & all(!abbr_is)) {
    cli_abort(c("Compound must available in open TG-GATEs database.",
                "x" = "{style_bold(col_red(backtick(comp_name)))} compound{?s} name/abbrerviation \\
                {?is/are} not available in open TG-GATEs database.",
                "i" = "Please find available compound name/abbrerviation by calling the function `compound()`.")
              , call = error_call)
  }
}

# compound name to abbreviation
comp_abbr <- function(comp_name, error_call = caller_env()) {
  if (is_empty(comp_name)) {
    cli_abort(c("{.var comp_name} must be non-empty.",
                "i" = "You have supplied an empty vector, please provide `name` or, `abbreviation` instred.")
              , call = error_call)
  }
  comp_tg <- compounds()
  comp_is <- comp_name %in% comp_tg$COMPOUND_NAME
  if (all(comp_is)) {
    idx <- match(comp_name, comp_tg$COMPOUND_NAME)
    abbr <- comp_tg$COMPOUND_ABBREVIATION[idx]
  }
  abbr_is <- comp_name %in% comp_tg$COMPOUND_ABBREVIATION
  if (all(abbr_is)) {
    abbr <- comp_name
  }
  if (any(comp_is) & any(abbr_is)) {
    cli_abort(c("Compound must be either in `name` or, `abbreviation`.",
                "x" = "Input {style_bold(col_blue(backtick(comp_name[comp_is])))} {?is/are} in compound name{?s} \\
                and {style_bold(col_green(backtick(comp_name [abbr_is])))} {?is/are} in compound abbreviation{?s}.",
                "i" = "Find available compound name/abbreviation by calling the function `compound()`.")
              , call = error_call)
  } else if (all(!comp_is) & all(!abbr_is)) {
    cli_abort(c("Compound must available in open TG-GATEs database.",
                "x" = "{style_bold(col_red(backtick(comp_name)))} compound{?s} name/abbrerviation \\
                {?is/are} not available in open TG-GATEs database.",
                "i" = "Please find available compound name/abbrerviation by calling the function `compound()`.")
              , call = error_call)
  }else if (any(!comp_is) & all(!abbr_is)) {
    if (any(!comp_is)) {
      if (sum(comp_is) > length(comp_name)/2) {
        idx <- match(comp_name[comp_is], comp_tg$COMPOUND_NAME)
        abbr <- comp_tg$COMPOUND_ABBREVIATION[idx]
        cli_alert_warning(c("Input {style_bold(col_red(backtick(comp_name[!comp_is])))} name{?s} of compound ",
                           "{?is/are} removed from the compound list, ",
                           "because {?this/these} input compound{?s} {?is/are} not available ",
                           "in open TG-GATEs database."),
                          wrap = TRUE)
      } else {
        cli_abort(c("Compound must available in open TG-GATEs database.",
                    "x" = "Input ({style_bold(col_red(backtick(comp_name[!comp_is])))}) name{?s} \\
                    of compound {?is/are} not available in open TG-GATEs database.",
                    "i" = "Please choose compound from the open TG-GATEs database,  \\
                    you can find compound list by calling the function `compound()`.")
                  , call = error_call)
      }
    }
  } else if (all(!comp_is) & any(!abbr_is)) {
    if (any(!abbr_is)) {
      if (sum(abbr_is) > length(comp_name)/2) {
        idx <- match(comp_name[abbr_is], comp_tg$COMPOUND_ABBREVIATION)
        abbr <- comp_tg$COMPOUND_ABBREVIATION[idx]
        cli_alert_warning(c("Input {style_bold(col_red(backtick(comp_name[!abbr_is])))} abbreviation{?s} of ",
                           "compound{?s} {?is/are} removed from the compound list, ",
                           "because {?this/these} compound{?s} {?is/are} not available ",
                           "in open TG-GATEs database."),
                          wrap = TRUE)
      } else {
        cli_abort(c("Compound must available in open TG-GATEs database.",
                    "x" = "Input ({style_bold(col_red(backtick(comp_name[!abbr_is])))}) \\
                    abbreviation{?s} of compound {?is/are} not available in open TG-GATEs database.",
                    "i" = "Please choose compound from the open TG-GATEs database,  \\
                    you can find available compound list by calling the function `compound()`.")
                  , call = error_call)
      }
    }

  }
  return(abbr)
}


# abbreviation or combination of abbreviation and compound name to compound name
abbr2name <- function(comp_abbr, error_call = caller_env()) {
  if (is_empty(comp_abbr)) {
    cli_abort(c("{.var comp_abbr} must be non-empty.",
                "i" = "You have supplied an empty vector, please provide `name` instred.")
              , call = error_call)
  }
  comp_tg <- compounds()
  comp_out <- comp_abbr %in% c(comp_tg$COMPOUND_ABBREVIATION, comp_tg$COMPOUND_NAME)
  if (any(!comp_out)) {
    cli_abort(c("Compound name must available in open TG-GATEs database.",
                "x" = "{style_bold(col_red(backtick(comp_abbr[!comp_out])))} compound{?s} name \\
                {?is/are} not available in open TG-GATEs database.",
                "i" = "Please find available compound name by calling the function `compound()`.")
              , call = error_call)
  }
  comp_is <- comp_abbr %in% comp_tg$COMPOUND_ABBREVIATION
  if (all(comp_is)) {
    idx <- match(comp_abbr, comp_tg$COMPOUND_ABBREVIATION)
    name <- comp_tg$COMPOUND_NAME[idx]
  }
  name_is <- comp_abbr %in% comp_tg$COMPOUND_NAME
  if (all(name_is)) {
    cli_alert_warning(c("Input {.var comp_abbr} {style_bold(col_red(backtick(comp_abbr)))} ",
                        "{?is/are} already in name{?s},  no need to convert."),
                      wrap = TRUE)
    name <- comp_abbr
  }
  if (any(name_is)) {
    cli_alert_warning(c("Input {.var comp_abbr} {style_bold(col_blue(backtick(comp_abbr[name_is])))} ",
                        "{?is/are} already in name{?s}"),
                      wrap = TRUE)
    rem_abbr <- comp_abbr[!name_is]
    rem_idx <- match(rem_abbr, comp_tg$COMPOUND_ABBREVIATION)
    rem_name <- comp_tg$COMPOUND_NAME[rem_idx]
    name <- comp_abbr
    name[!name_is] <- rem_name
  }
  return(name)
}


# Find data type from attribute data
get_dtype <- function(attr_df, error_call = caller_env()) {
  species <- unique(attr_df$species)
  test_type <- unique(attr_df$test_type)
  organ <- unique(attr_df$organ_id)
  sin_rep <- unique(attr_df$sin_rep_type)
  if (length(species) > 1) {
    cli_abort(c("More than one `species` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} species {?is/are}  \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `species` to get appropiate data type.")
              , call = error_call)
  }
  if (length(test_type) > 1) {
    cli_abort(c("More than one `test_type` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} test type{?s} {?is/are}  \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `test_type` to get appropiate data type.")
              , call = error_call)
  }
  if (length(organ) > 1) {
    cli_abort(c("More than one `organ` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} organ{?s} {?is/are}  \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `organ` to get appropiate data type.")
              , call = error_call)
  }
  if (length(sin_rep) > 1) {
    cli_abort(c("More than one `sin_rep_type` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} experiment type{?s} {?is/are} \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `sin_rep_type` to get appropiate data type.")
              , call = error_call)
  }
  dtype <- case_when(
    species == "Rat" & test_type == "in vivo" & organ == "Liver" & sin_rep == "Single" ~ "R1LS",
    species == "Rat" & test_type == "in vivo" & organ == "Liver" & sin_rep == "Repeat" ~ "R1LR",
    species == "Rat" & test_type == "in vivo" & organ == "Kidney" & sin_rep == "Single" ~ "R1KS",
    species == "Rat" & test_type == "in vivo" & organ == "Kidney" & sin_rep == "Repeat" ~ "R1KR",
    species == "Rat" & test_type == "in vitro" & organ == "Liver" & is.na(sin_rep) ~ "R2L",
    species == "Human" & test_type == "in vitro" & organ == "Liver" & is.na(sin_rep) ~ "H2L"
  )
  return(dtype)
}


select_pcol <- function(pmat, lab = NULL, error_call = caller_env()) {
    if (is.null(lab)) {
      test_pmat <- pmat
    } else {
      lab <- gsub(" ", "", lab)
      all_lab <- c("group",  "compound", "dose", "time")
      unlist_lab <- unlist(strsplit(lab, "[^(A-Za-z)]"))
      if (any(unlist_lab %in% "")) {
        unlist_lab <- unlist_lab[-which(unlist_lab == "")]
      }
      if (!all(unlist_lab %in% all_lab)) {
        miss_lab <- unlist_lab[!unlist_lab %in% all_lab]
        remin_lab <- all_lab[!all_lab %in% unlist_lab]
        cli_abort(c("The name of stage used in test level must be a valid stage.",
                    "x" = "Input {style_bold(col_red(backtick(miss_lab)))} {?is/are} not a valid stage.",
                    "i" = "Please use {style_bold(col_red(backtick(remin_lab)))} instread."),
                  call = error_call)
      }
      join_sym <- unlist(strsplit(gsub("[^[:punct:]S]", "", lab), ""))
      if (!all(join_sym %in% "&")) {
        wrong_sym <- join_sym[!join_sym %in% "&"]
        cli_abort(c("Test level condition contain wrong symbol.",
                    "x" = "Input {style_bold(col_red(backtick(wrong_sym)))} {?is/are} not a valid symbol.",
                    "i" = "Please use `&` instread."),
                  call = error_call)
      }
      test_pmat <- as.matrix(pmat[, unlist_lab])
    }
    return(test_pmat)
  }
# pmat <- matrix(rnorm(20),ncol = 4)
# colnames(pmat) <- c("group",  "compound", "dose", "time")
# rownames(pmat) <- LETTERS[1:5]
# select_pcol(pmat, lab="group")


# up-doen regulation calculator

up_down <- function(...,
                    probes,
                    expr_data,
                    attr_data,
                    multicore = FALSE,
                    store = FALSE,
                    output_dir = missing_arg(),
                    error_call = caller_env()) {
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
  lev_str <- expr_str(comps_gr, attr_data = comp_dict)
  design_mat <- design_matrix(lev_str)
  if (is.null(names(comps_gr))) {
    names(comps_gr) <- paste("Group", LETTERS[1:length(comps_gr)], sep = "_")
  }
  avg_gfc <- expr_data %**% design_mat$qAlfa
  colnames(avg_gfc) <- names(comps_gr)
  rownames(avg_gfc) <- rownames(expr_data)
  up_downrg <- apply(avg_gfc, 1, function(x) ifelse(x[1]>x[2], "Upregulated", "Downregulated"))
  updown_df <- tibble::as_tibble(avg_gfc, rownames = "probe_id") %>%
    dplyr::mutate(Expression = factor(up_downrg,levels = c("Upregulated", "Downregulated")))
  return(updown_df)
}


split_two <- function(sentence) {
  words <- strsplit(sentence, "\\s+")[[1]]
  num_words <- length(words)

  if (num_words < 2) {
    first_half <- sentence
    second_half <- ""
  } else {
    mid_point <- ceiling(num_words / 2)
    first_half <- paste(words[1:mid_point], collapse = " ")
    second_half <- paste(words[(mid_point + 1):num_words], collapse = " ")
  }
  return(list(first_half = first_half, second_half = second_half))
}


split_multiple <- function(sentence, nparts) {
  words <- strsplit(sentence, "\\s+")[[1]]
  num_words <- length(words)

  if (num_words < nparts) {
    cli_abort(glue::glue("The sentence should contain at least {nparts} words."))
  }

  words_per_part <- ceiling(num_words / nparts)

  sentence_parts <- vector("list", nparts)

  for (i in 1:nparts) {
    start_idx <- (i - 1) * words_per_part + 1
    end_idx <- min(start_idx + words_per_part - 1, num_words)
    sentence_parts[[i]] <- paste(words[start_idx:end_idx], collapse = " ")
  }

  return(sentence_parts)
}


split2_merge <- function(sentences) {
  split_x <- lapply(sentences, function(x) split_two(x))
  max_char <- max(sapply(unlist(split_x), nchar))
  merge_split <- unlist(lapply(split_x,
                        function(x) {ifelse(x[[2]] != "",
                                            paste(unlist(x),collapse = "\n"), x[[1]])}))
  merge_split[nchar(sentences) < max_char] <- sentences[nchar(sentences) < max_char]
  return(merge_split)
}

title_tag <- function(title = NULL, title_size = rel(2), note = NULL,
                      note_size = rel(0.5), tag = NULL, tag_size = rel(1)) {
  if (!is.null(tag)) {
    text(x = -1,
         y = 1 ,
         labels = tag,
         cex = tag_size,
         font = 2)}
  if (!is.null(title)) {
    title(main = title,
          cex.main = title_size,
          font.main = 2
    )
  }
  if (!is.null(note)) {
    text(x = 1,
         y = -1 ,
         labels = note,
         pos = 2,
         cex = note_size
    ) #, family= "Times New Roman")
  }
}

handl_group <- function(..., error_call = caller_env()) {
  comps_gr <- list2(...)
  if (length(comps_gr) == 1 && is_bare_list(comps_gr[[1]])) {
    comps_gr <- comps_gr[[1]]
  }
  arg_checker <- unlist(lapply(comps_gr, function(x) inherits(x, c("character", "factor", "list"))))
  right_arg <- comps_gr[arg_checker]
  right_class <- unlist(lapply(right_arg, class))
  if (length(unique(right_class))>1) {
    right_name <- sapply(substitute(list(...))[-1], as.character)[arg_checker]
    cli_abort(c("Invalid inpute of compound groups.",
                "x" = "You have supplied the object{?s} {paste(style_bold(col_red(right_name)), style_bold(col_blue(right_class)), sep = ' as a ')}
                class{?, respectively}.",
                "i" = "Please provide the compound groups as either a list or multiple vectors containing compound names or abbreviations."),
              wrap = TRUE, call = error_call)
  }
  # wrong_arg <- which(!arg_checker)
  wrong_arg <- sapply(substitute(list(...))[-1], as.character)[!arg_checker]
  if (!all(arg_checker)) {
    arg_nm <- paste("argument_name", LETTERS[1:length(wrong_arg)])
    cli_abort(c("Invalid use of the argument in the function.",
                "x" = "You have supplied object{?s} {style_bold(col_red(wrong_arg))}  without assigning
                {?it/them} as an argument{?s}. However, {?it/they} {?does/may} not contain any compound names or abbreviations.",
                "i" = "Please provide the argument explicitly by supplying the compound names or abbreviations correctly.
                Alternatively, if applicable, you can assign the appropiate argument{?s} in the function as:
                {paste(style_bold(col_green('argument_name')), style_bold(col_red(wrong_arg)), sep = ' = ')}."),
              wrap = TRUE, call = error_call)
  }
  if (is.null(names(comps_gr))) {
    names(comps_gr) <- paste("Group", LETTERS[1:length(comps_gr)], sep = "_")
  }
  return(comps_gr)
}


rm_msg <- function(cat_msg) {
  cat("\r")
  cat(paste0(rep(" ", nchar(cat_msg)), collapse = ""))
  cat("\r")
}



get_class <- function(comps_gr, compounds) {
  class_assignment <- sapply(compounds, function(element) {
    matched_clusters <- sapply(comps_gr, function(cluster) element %in% cluster)
    match(TRUE, matched_clusters)
  })

  return(class_assignment)
}



round_up <- function(x, digits = 2) {
  multiplier <- 10^digits
  ifelse(x > 0, ceiling(x * multiplier) / multiplier, floor(x * multiplier) / multiplier)
}


get_organism <- function(attr_data) {
  chip <- unique(attr_data$arr_design)
  if(identical(chip, "Rat230_2")) {
    organism = "rat"
  }else if(identical(chip, "HG-U133_Plus_2")) {
    organism = "human"
  }
  return(organism)
}

get_genename <- function(probes, organism) {
  genes <- probes2genes(probes, organism) %>%
    mutate(SYMBOL = dplyr::coalesce(SYMBOL, PROBEID))
  gene_names <-  genes$SYMBOL
  # sample_expr <- sample_expr[!duplicated(rownames(sample_expr)), ]
  return(gene_names)
}
