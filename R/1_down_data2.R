#' List of compounds available in TG-GATES database
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The `compounds()` function is used to get a list of available compounds and their abbreviations
#' that will use to download the raw gene expression data.
#'
#' The `compounds()` function is used to get a list of available compounds in
#' [open TG-GATEs](https://toxico.nibiohn.go.jp/english/) database. The compound name or its abbreviation is
#' used to download and processing affymetric CEL image of samples (combinations of `compound`, `dose`,
#' `timepoint` and `replications`) and their attributes data.
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return
#' A data frame of 2 columns
#'   * Column-1: Compound name
#'   * Column-2: Compound abbreviations
#'   * Column-3: Compound number
#' @export
#'
#' @examples
#' compounds()
load(file="E:/TGGATES_v0.0/data/comp_data.rda") # delete it when package ready
compounds <- function() {
  comp_data %>%
    dplyr::select(c("COMPOUND_NAME", "COMPOUND_ABBREVIATION", "COMPOUND_NO"))
}


#' Download link of a compound from open TG-GATEs
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `compound_link()` is used to get the download link of query compound.
#'
#' @param comp_name Name of the compound or its abbreviation.
#' @param data_type Compound data type. Compound was administrated to **species (Rat/Human)**
#'    following a experimental condition **test type (in-vivo/in-vitro)** and then
#'    extracted data from **organ (Liver/Kidney)** for different **dose type (Singe/Repeated)**.
#'    Therefore, the `data_type` follow  the sequence `species -> test type -> organ -> dose type`.
#'    Available data types:
#'
#'    * `"R1LS"`: For Rat in-vivo liver Single data.
#'    * `"R1LR"`: For Rat in0-vivo liver Repeat data.
#'    * `"R1KS"`: For Rat in-vivo kidney Single data.
#'    * `"R1KR"`: For Rat in-vivo kidney Repeat data.
#'    * `"R2L"`: For Rat in-vitro data.
#'    * `"H2L"`: For Human in-vitro data.
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return
#' A download link
#' @export
#'
#' @examples
#' compound_link("2NF", "R1LS")
compound_link <- function(comp_name, data_type = c("R1LS", "R1LR", "R1KS", "R1KR", "R2L", "H2L")) {
  # data_type <- match.arg(data_type)
  #formal_args <- formals(sys.function(sys_par <- sys.parent()))
  #data_types <- eval(formal_args[["data_type"]], envir = sys.frame(sys_par))
  data_types <- c("R1LS", "R1LR", "R1KS", "R1KR", "R2L", "H2L")
  if (!data_type %in% data_types) {
    cli_abort(c(
      "Wrong data type",
      "x" = "The {.arg data_type} {style_bold(col_red(backtick(data_type)))} you provided is not recognized",
      "i" = "Please chosse either {style_bold(col_green(backtick(data_types)))} instead"
    ))
  }
  comp_name <- comp_abbr(comp_name)
  col_index <- switch(data_type,
                      R1LS = match("Rat...in.vivo...Liver...Single", names(comp_data)),
                      R1LR = match("Rat...in.vivo...Liver...Repeat", names(comp_data)),
                      R1KS = match("Rat...in.vivo...Kidney...Single", names(comp_data)),
                      R1KR = match("Rat...in.vivo...Kidney...Repeat", names(comp_data)),
                      H2L  = match("Human...in.vitro", names(comp_data)),
                      R2L  = match("Rat...in.vitro", names(comp_data)))

  # if (identical(data_type, "R1LS")) {
  #   data_type <- "Rat...in.vivo...Liver...Single"
  #   col_index <- match(data_type, names(comp_data))
  # } else if (identical(data_type, "R1LR")) {
  #   data_type <- "Rat...in.vivo...Liver...Repeat"
  #   col_index <- match(data_type, names(comp_data))
  # } else if (identical(data_type, "R1KS")) {
  #   data_type <- "Rat...in.vivo...Kidney...Single"
  #   col_index <- match(data_type, names(comp_data))
  # } else if (identical(data_type, "R1KR")) {
  #   data_type <- "Rat...in.vivo...Kidney...Repeat"
  #   col_index <- match(data_type, names(comp_data))
  # } else if (identical(data_type, "H2L")) {
  #   data_type <- "Human...in.vitro"
  #   col_index <- match(data_type, names(comp_data))
  # } else if (identical(data_type, "R2L")) {
  #   data_type <- "Rat...in.vitro"
  #   col_index <- match(data_type, names(comp_data))
  # } else {
  #   cli_abort(c(
  #     "Wrong data type",
  #     "x" = "The {.var data_type} ({data_type}) you provided is not recognized",
  #     "i" = "Please chosse either `R1LS`, `R1LR`, `R1KS`, `R1KR`, `H2L` or `R2L` instead"
  #   )
  #   )
  # }
  row_index <- match(comp_name, comp_data$COMPOUND_ABBREVIATION)
  drug_link <- comp_data[row_index, col_index]
  if (identical(drug_link,"")) {
    other_dtype <- data_types[!data_types %in% data_type]
    cli_abort(c(
      "No link found",
      "x" = "There is no data for compound {style_bold(col_red(backtick(comp_name)))} in {style_bold(col_br_red(backtick(data_type)))}",
      "i" = "Please chosse {style_bold(col_green(backtick(other_dtype)))} instead"
    ))
  }
  return(drug_link)
}

#' Download CEL image files and attribute data of compound from open TG-GATEs
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `download_rawdata()` is used to download affymetric CEL image files and the attribute
#' (which contain data dictionary and phenotypic information) data of a compound. A folder will create in
#' your output directory with name `compound_name.data_type`, a sub-folder named `celfiles` is create
#' under this parent folder that contain all the CEL image files. However, a `.tsv` file of attributes
#' information is also create to this parent folder.
#'
#' The `download_rawdata()` is used to download affymetric CEL image files of compound from
#' [open TG-GATEs](https://toxico.nibiohn.go.jp/english/) database. All the data will save in your selected
#' directory. Since the expression of genes is measured at different doses and time points or combination thereof,
#' resulting a set of `.CEL` files that will download to a folder named `celfiles` under the parent folder
#' `compound_name.data_type`(compound name dot data type). A attribute `.tsv` file contain all phenotypic
#' (including biochemical and biochemistry) and data dictionary information also saved in parent directory.
#' Keep in mind that the download speed is depends on your internet speed.
#'
#' @inheritParams compound_link
#' @param output_dir Output directory where you want to download raw data. Current working directory will
#'    be use if output directory is missing.
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A folder with `.CEL` files and a `.tsv` attribute file
#' @export
#'
#' @examples
#' \donttest{
#' download_rawdata(comp_name="2NF", data_type="R1LS")
#' }
download_rawdata <- function(comp_name, data_type, output_dir = missing_arg()) {
  # change the download need for `download.file` function timeout since the files are big
  # opts <- options()
  # options(timeout=1000)
  # on.exit(options(opts))
  output_dir <- destination(output_dir)
  if (!curl::has_internet()) {
    cli_abort(c("No internet connection.",
                "i" = "Please connect an internet and try again."))
  }
  comp_link <- compound_link(comp_name = comp_name, data_type = data_type)
  if (!RCurl::url.exists(comp_link)) {
    cli_abort(c("Open TG-GATEs FTP server not responding.",
                "i" = "Please try again later."))
  }
  temp_comp <- gsub("^.*/", "", comp_link)
  sp_comp <- unlist(strsplit(temp_comp, "\\."))
  dt_type <- glue_collapse(sp_comp[-c(1, length(sp_comp))],"-")
  out_path <- glue::glue("{output_dir}/{temp_comp}")
  start_time <- Sys.time()
  #utils::download.file(comp_link, out_path, quiet = TRUE)
  #curl::curl_download(comp_link, out_path)
  cli_alert_info("Downloading {style_bold(col_red(backtick(sp_comp[1])))} ({col_br_red(dt_type)}) data:")
  #httr::GET(comp_link, httr::write_disk(out_path, overwrite=TRUE), httr::progress("down"))
  #utils::download.file(url = comp_link, destfile = out_path, quiet = TRUE)
  httr::GET(comp_link,httr::write_disk(out_path, overwrite = TRUE), httr::progress())
  end_time <- Sys.time()
  req_time <- difftime(end_time, start_time, units = "secs")[[1]]
  file_size <- file.info(out_path)$size
  zip_file <- list.files(path = output_dir, pattern = temp_comp, full.names = TRUE)
  utils::unzip(zip_file, exdir = output_dir)
  rm_zip <- file.remove(zip_file)
  #sp_comp <- unlist(strsplit(temp_comp, "\\."))
  #dt_type <- glue_collapse(sp_comp[-c(1, length(sp_comp))],"-")
  if(rm_zip) {
    cli_alert_success(c("Downloaded {prettyunits::pretty_bytes(file_size)}",
                             " in {prettyunits::pretty_sec(as.numeric(req_time))}"))
  }else {
    cli_alert_success("{style_bold(col_red(backtick(sp_comp[1])))} ({col_br_red(dt_type)}) data download has unsuccessful")
  }
  path <- gsub(".zip$","" , out_path)
  return(path)
}

#' CEL images to gene expression matrix
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `cel2exprs()` is used to extract gene expression data from CEL files. The
#' dimension of output expression matrix depends on the number of CEL files and the probes in CDF.
#'
#' `cel2exprs()` is use core functions [`mas5`][affy::mas5()] and [`rma`][affy::rma()] from `affy` package
#' to quantify mRNA intensity of probes/genes. Note that expression measure is given in log base2 scale.
#'
#' @param cel_files A vector of full path of CEL files.
#' @param organism Sample organism used for experiment either "human" or "rat".
#' @param method Method used for normalization
#' @param ... Arguments used for [`mas5`][affy::mas5()] or [rma][affy::rma()], e.g., `normalize`, `background` etc.
#'
#' @family Open TG-GATEs data download helpers
#'
#' @seealso
#'    [affy::mas5()] for Affymetrix version 5 (MAS5) model,
#'    [affy::rma()] for Robust Means Analysis (RMA) model
#'
#' @return
#' A matrix of raw gene expression data.
#'   * Probes in row
#'   * CEL ID in column
#' @export
#'
#' @examples
#' \donttest{
#' file_path <- download_rawdata("2NF", "R1LS")
#' cel_file <- paste(file_path, "celfiles", sep = "/")
#' cel_files <- list.files(cel_file, full.names = TRUE)
#' expr_data <- cel2exprs(cel_files)
#' head(expr_data)
#' }
cel2exprs <- function(cel_files, organism = "rat", method = "rma", ...)
  {
  if (identical(organism, "human")) {
    CDF = "hgu133plus2cdf"
  } else if(identical(organism, "rat")) {
    CDF = "rat2302cdf"
  } else {
    cli_abort(c(
      "Wrong organism name",
      "x" = "The name of {.var organism} ({style_bold(col_blue(backtick(organism)))}) you provided is not recognized",
      "i" = "Please chosse either `human` or `rat` instead"
      ))
  }
  files <- tools::file_ext(cel_files)
  cel_path <- cel_files[which(files == "CEL")]
  if (!all(files == "CEL")) {
    cli_warn(c(
      "Some of the path not contain any CEL file:",
      i = "The following path should be removed",
      format_bullets_raw(rlang::set_names(cel_files[!files == "CEL"], "x"))
      ))
  }
  read_cel <- affy::ReadAffy(filenames = cel_path, cdfname = CDF)
  if (identical(method, "mas5")) {
    expr_meas <- affy::mas5(read_cel, verbose = FALSE, ...)
    expr <- Biobase::exprs(expr_meas)
    expr <- log2(expr)
  } else if (identical(method, "rma")) {
    expr_meas <- affy::rma(read_cel, verbose = FALSE, ...)
    expr <- Biobase::exprs(expr_meas)
  }
  return(expr)
}

#' Process gene expression and meta data
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `process_celfile()` is used to create gene expression and metadata of samples
#' from CEL  image and attribute files respectively.
#'
#'
#' @param comp_path Full path of compound that contain a folder `celfile` and a `.tsv` attribute file
#' @param log2FC A boolean.
#'   * `TRUE` (the default): Output gene expression data will be Log2 fold-change.
#'   * `FALSE`: Raw gene expression data.
#' @param multicore Use for parallel computing? If `TRUE` (the default), Computation will use parallel
#'   processing. If `FALSE`, no palatalization will use. To apply parallel computing use either `TRUE`,
#'   or an integer to specify number of workers. For `TRUE`, computation will left 1 core for other uses.
#'   Its better to leave one core for system use. For example, if your machine has 8 core use 7 for your computation.
#' @param output_dir A path where you want to download the data. If path is missing, current working
#'   directory will be used.
#' @inheritParams cel2exprs
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A list.
#'   * A matrix of gene expression measure, probes in row and sample (CEL-ID) in column.
#'   * A data frame of metadata about sample and phenotypic information .
#' @export
#'
#' @examples
#' \donttest{
#' file_path <- download_rawdata("2NF", "R1LS")
#' proc_expr <- process_celfile(file_path)
#' head(proc_expr$gexpr_data)
#' head(proc_expr$meta_data)
#' }
process_celfile <- function(comp_path,
                            fc = TRUE,
                            method = "rma",
                            multicore = FALSE,
                            output_dir = missing_arg(),
                            ...) {
  output_dir <- destination(output_dir)
  start_time <- Sys.time()
  sp_path <- unlist(strsplit(comp_path, "\\."))
  sp_comp <- unlist(strsplit(sp_path[1], "/"))
  comp_nm <- sp_comp[length(sp_comp)]
  comp_typ <- glue_collapse(sp_path[-1],"-")
  start_msg <- "Processing {style_bold(col_red(backtick(comp_nm)))} ({col_br_red(comp_typ)}) data........."
  cat(glue(start_msg))
  attr_file <- paste(comp_path, "Attribute.tsv", sep = "/")
  #comp_dict <- readr::read_tsv(file = attr_file, show_col_types = FALSE)  # not work for rat data
  comp_dict <- utils::read.table(file = attr_file, sep = '\t', header = TRUE) # not work for human data
  comp_dict <- comp_dict %>%
    dplyr::rename_with(~ tolower(gsub(" ", "_", .x, fixed = TRUE))) %>%
    dplyr::rename(
      compound_abbr = "compound.abbr.",
      time_level = "sacri_period",
      sample_id = "individual_id"
    ) %>%
    dplyr::filter(.data$barcode != "No ChipData") %>%
    dplyr::mutate_at(c("compound_name", "compound_abbr", "dose_level", "time_level"), as.factor) %>%
    dplyr::relocate(
      c("compound_name", "compound_abbr", "dose_level", "time_level",
        "sample_id", "arr_design", "species", "test_type", "organ_id", "sin_rep_type"),
      .after = "barcode")
  chip <- unique(comp_dict$arr_design)
  if(identical(chip, "Rat230_2")) {
    organism = "rat"
  }else if(identical(chip, "HG-U133_Plus_2")) {
    organism = "human"
  }
  cel_file <- paste(comp_path, "celfiles", sep = "/")
  file_path <- list.files(path = cel_file, pattern = "*.CEL", full.names = TRUE)
  if(is.logical(multicore)) {
    if(multicore ) {
      multicore <- start_parallel(multicore)
      stop_cluster <- TRUE
    } else {
      multicore <- stop_cluster <- FALSE
    }
  }else {
    stop_cluster <- if(inherits(multicore, "cluster")) FALSE else TRUE
    multicore <- start_parallel(multicore)
  }
  on.exit(if(multicore & stop_cluster)
    stop_parallel(attr(multicore, "cluster")))
  `%DO%` <- if(multicore) foreach::`%dopar%` else foreach::`%do%`
  if(!multicore) {
    # cli_alert_warning(c("Processing without parallel computing will take much time. ",
    #                   "Please use `multicore = TRUE` to parallelize your task."))
    expr_data <- cel2exprs(file_path, organism = organism, method = method, ...)
  }else {
    ii <- attr(multicore, "cores")
    # Set condition if ii>length(file_path), then use length(file_path) cores only.
    par_idx <- parallel::splitIndices(length(file_path), ii)
    fun_call <- as.call(c(list(quote(foreach::foreach), i = seq_len(ii)), .combine ="cbind"))
    .fun <- eval(fun_call)
    #cli_alert_warning(c("Processing will takes time."))
    expr_data <- .fun %DO% cel2exprs(file_path[par_idx[[i]]], organism = organism, method = method, ...)
  }
  colnames(expr_data) <- gsub(".CEL", "", colnames(expr_data))
  timepnt <- levels(comp_dict$time_level)
  temp_dict <- comp_dict %>%
    dplyr::filter(.data$dose_level != "Control")
  if(fc) {
    temp_expr <- list()
    for(i in 1:length(timepnt)) {
      const_dict <- comp_dict %>%
        dplyr::filter(.data$dose_level == "Control" & .data$time_level == timepnt[i])
      treat_dict <- comp_dict %>%
        dplyr::filter(.data$dose_level != "Control" & .data$time_level == timepnt[i])
      const_expr <- expr_data[, colnames(expr_data) %in% const_dict$barcode]
      treat_expr <- expr_data[, colnames(expr_data) %in% treat_dict$barcode]
      const_med <- apply(const_expr, 1, stats::median)
      temp_expr[[i]] <- treat_expr - const_med
    }
    raw_expr <- do.call(cbind, temp_expr)
    match_col <-  match(temp_dict$barcode , colnames(raw_expr))
    match_col <- match_col[!is.na(match_col)]
    exprs_data <-raw_expr[, match_col]
    temp_dict <- temp_dict %>%
      mutate(fc = TRUE, .before = "arr_design")
  }else {
    match_col <-  match(temp_dict$barcode, colnames(expr_data))
    match_col <- match_col[!is.na(match_col)]
    exprs_data <- expr_data[, match_col]
  }
  out_data <- list(exprs_data, temp_dict)
  names(out_data) <- c("gexpr_data", "meta_data")
  end_time <- Sys.time()
  if(!multicore) {
    req_time <- difftime(end_time, start_time, units = "secs")[[1]]
    rm_msg(start_msg)
    cli_alert_success(c("Processing of {style_bold(col_red(backtick(comp_nm)))}",
                    "({col_br_red(comp_typ)}) data has been completed in ",
                    "{prettyunits::pretty_sec(as.numeric(req_time))}."))
  }else {
    req_time <- difftime(end_time, start_time, units = "secs")[[1]]
    rm_msg(start_msg)
    cli_alert_success(c("Processing of {style_bold(col_red(backtick(comp_nm)))}",
                        "({col_br_red(comp_typ)}) data has been completed ",
                        "in {prettyunits::pretty_sec(as.numeric(req_time))} using {ii} worker{?s}."))
  }
  return(out_data)
}

#' Process a set of compounds data from open TG-GATEs
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' `process_compaunds()` process a set of compounds from open TG-GATEs database. It process all CEL image and
#' attribute files of compound and generate two CSV file `expression.csv` and `metadata.csv`, respectively. if
#' you set the flag `store = TRUE`, both `.csv` files will save in your output directory.
#'
#' @param files A vector of compounds path which contain a folder `celfile` and a `.tsv` attribute file.
#' @inheritParams process_celfile
#' @param store Save resulting data in the output directory? IF `TRUE` (the default), two `.csv`files
#'  of gene expression and metadata data will save in the output directory.
#'
#' @family Open TG-GATEs data download helpers
#' @return A list.
#'   * A matrix of gene expression measure of all compounds, probes in row and sample (CEL-ID) in column.
#'   * A data frame of attribute about sample and phenotypic information of compounds.
#' @export
#'
#' @examples
#' \donttest{
#' files <- lapply(list("2NF", "NMOR"),download_rawdata, data_type = "R1LS")
#' tgx_data <- suppressMessages(process_compaunds(files))
#' head(tgx_data$expression)
#' head(tgx_data$attribute)
#' }
process_compaunds <- function(files,
                              fc = TRUE,
                              method = "rma",
                              multicore = FALSE,
                              store = FALSE,
                              output_dir = missing_arg(),
                              ...) {
  files <- as.list(files)
  exprs_list <- lapply(files, process_celfile, fc = fc,
                       multicore = multicore, method = method, output_dir = output_dir, ...)
  expr_mat <- do.call(cbind, lapply(exprs_list, "[[",1))
  dict_df <- do.call(rbind, lapply(exprs_list, "[[",2))
  if (store) {
    output_dir <- destination(output_dir)
    expr_mat %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    dict_df %>%
      readr::write_csv(file = paste(output_dir, "metadata.csv", sep = "/"))
  }
  tgx_data <- list(expr_mat, dict_df)
  names(tgx_data) <- c("expression", "metadata")
  return(tgx_data)
}

#' Download and process open TG-GATES data
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' `get_TGGATEs()` directly download and process compounds data from the
#' [open TG-GATEs](https://toxico.nibiohn.go.jp/english/) database.
#'
#' The function `get_TGGATEs()` is used to download and process compounds raw data. This function
#' first download raw CEL images data of compounds and then process them into expression measure.
#' All the compounds will download to the output directory with the folder name
#' `compound_name.datatype` (compound name dot data type). The raw CEL images are stored in a
#' sub-folder `celfiles` under the parent folder `compound_name.datatype`. However, a attribute
#' TSV file also download for each compounds. Finally, two CSV file `expression.csv` and
#' `metadata.csv` for gene expression and attribute data, respectively are saved in output directory.
#' A list of available compounds can be found by calling the function [compounds()].
#'
#' @param compds A vector of compounds name or its abbreviation
#' @inheritParams compound_link
#' @inheritParams process_compaunds
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A list.
#'   * A matrix of gene expression measure of all compounds, probes in row and sample (CEL-ID) in column.
#'   * A data frame of attribute about sample information of compounds.
#' @export
#'
#' @examples
#' \donttest{
#' compds <- c("2NF", "NMOR")
#' otg_data <- get_TGGATEs("2NF", data_type="R1LS",output_dir="D:/project1_TGx/data_processing4")
#' head(otg_data$expression)
#' head(otg_data$metadata)
#' }

get_TGGATEs <- function(compds,
                        data_type,
                        fc = TRUE,
                        method = "rma",
                        multicore = TRUE,
                        store = FALSE,
                        output_dir = missing_arg(),
                        ...) {
  compds <- as.list(compds)
  comp_name <- lapply(compds, comp_abbr)
  raw_files <- lapply(
    compds,
    download_rawdata,
    data_type = data_type,
    output_dir = output_dir)
  tg_data <- process_compaunds(
    raw_files,
    fc = fc,
    multicore= multicore,
    method = method,
    output_dir = output_dir,
    store = store,
    ...)
  return(tg_data)
}
