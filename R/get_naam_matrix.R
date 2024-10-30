

#' Structure of experiment to test a toxic endpoint
#'
#' The function `expr_str()` is used to find the number of levels in dose, time point and replication
#' of a list of compounds in different groups. `expr_str()` uses the function `comp_str()` to find
#' the structure of levels in compound and concatenate them to create structural levels of experiment.
#'
#' @param comps_gr A list of compounds in group
#' @param attr_data Attribute data of compounds
#'
#' @return
#' A list total levels in different stage or a class of `tggates` object .
#'  * group
#'  * compound
#'  * dose
#'  * time
#'  * replication
#' @export
#'
#' @examples
#' \dontrun{
#' attr_data <- readr::read_csv("D:/project1_TGx/data_processing/data_dict.csv")
#' comps_gr <- list("APAP", "BBZ")
#' expr_str(comps_gr, attr_data = attr_data)
#' }
level_str2 <- function(chemicals, attr_data, error_call = caller_env()) {
  comps <- comp_abbr(chemicals, error_call = error_call)
  comp_lev <- length(chemicals)
  comps_list <- as.list(comps)
  zz <- lapply(comps_list, comp_str, attr_data = attr_data)
  dose_lev <- do.call(c,lapply(zz, "[[", 1))
  time_lev <- do.call(c,lapply(zz, "[[", 2))
  rep_lev <- do.call(c,lapply(zz, "[[", 3))
  comps_lev <- list(comp_lev, dose_lev, time_lev, rep_lev)
  names(comps_lev) <- c("compound", "dose", "time", "replication")
  structure(comps_lev, class = "NAAM")
  #return(comps_lev)
}


#' Design matrix for multistage hierarchical analysis of variance (MHANOVA) model
#'
#' The function `design_matrix()` is used to create design and other useful matrices for
#' MHANOVA model.
#'
#' @param lev_str structure of the experiment, a list of total number of levels in each
#'   stage or an `tggates` object of experimental structure. It can be obtained using
#'   the function `expr_str()`.
#'
#' @return a list of design matrices and others information in the experiment
#'   * `a`, total number of groups.
#'   * `b.`, total number of compound(s) in each group.
#'   * `c..`, total number of doses in each group and compound combination.
#'   * `d...`, total number of time points in each group, compounds and dose combination.
#'   * `N`, total number of samples.
#'   * `n_i...`, total number of samples in each group.
#'   * `n_ij..`, total number of samples in each compound within group.
#'   * `n_ijk.`, total number of samples in each dose within compound and group.
#'   * `n_ijkl`, total number of samples in each time-point within dose, compound and group.
#'   * `X1`, design matrix for group stage for MHANOVA model.
#'   * `X2`, design matrix for compound stage for MHANOVA model.
#'   * `X3`, design matrix for dose stage for MHANOVA model.
#'   * `X4`, design matrix for time stage for MHANOVA model
#'   * `QA`, information matrix used to calculate sum of square in group stage.
#'   * `QB`, information matrix used to calculate sum of square in compound stage.
#'   * `QC`, information matrix used to calculate sum of square in dose stage.
#'   * `QD`, information matrix used to calculate sum of square in time-point stage.
#'   * `R`, information matrix used to calculate sum of square error.
#'   * `QAlfa`, information matrix used to calculate design matrix `X1`.
#'   * `qAlfa`, information matrix used to calculate average expression at group level.
#'   * `QBeta`,
#'   * `qBeta`, information matrix used to calculate average expression at group level.
#'   * `QGama`,
#'   * `qGama`, information matrix used to calculate average expression at group level.
#'   * `QDelta`,
#'   * `qDelta`, information matrix used to calculate average expression at group level.
#'   * `Group`, group levels in sample
#'   * `Compound`, compound levels in sample.
#'   * `Dose`, dose levels in sample.
#'   * `Time`, time levels in sample.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' attr_data <- readr::read_csv("D:/project1_TGx/data_processing/data_dict.csv")
#' comps <- list("APAP", "BBZ")
#' lev_str <- expr_str(comps, attr_data = attr_data)
#' xx <- design_matrix(lev_str)
#' }
get_matrix <- function(lev_str) {
  if (!is.list(lev_str) & !inherits(lev_str, "tggates")) {
    cli_abort(c("{.var lev_str} must be object of class `tggates` or `list`.",
                "x" = "The class {style_bold(col_cyan(backtick(class(lev_str))))} of \\
                {.var lev_str} is not supported.",
                "i" = "Please provide {.var lev_str} as a list of size 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."))
  } else if (length(lev_str)  !=  4) {
    cli_abort(c("The length of {.var lev_str} must be 5.",
                "i" = "You have supplied a list {.var lev_str} of size {length(lev_str)}, please
                make sure {.var lev_str} has a length of 5
                (for `compound`, `dose`, `time` and `replication`)."))
  }
  b_i <- lev_str[[1]]
  c_ij <- lev_str[[2]]
  d_ijk <- lev_str[[3]]
  n_ijkl <- lev_str[[4]]
  N <-  sum(lev_str[[4]])
  Compound = Dose = Time = n_ijk. = n_ij.. = n_i... <- NULL
  x1 = x2 = x3 = Qa = Qb = Qc <- matrix(nrow=0, ncol=0)
  qBeta1 = qBeta2 = qGama1 = qGama2 = qDelta1 = qDelta2 <-  matrix(nrow=0, ncol=0)
    ni... <- 0
    for (j in 1 : b_i) {
      nij.. <- 0
      for (k in 1 : c_ij[j]) {
        nijk. <- 0
        for (l in 1 : d_ijk[k]) {
          nijk. <- nijk. + n_ijkl[l]
          x3 <- direct_sum(x3, matrix(rep(1, n_ijkl[l]), ncol = 1))
          qDelta1 <- direct_sum(qDelta1, matrix(rep(1, n_ijkl[l]), ncol = 1) / n_ijkl[l])
          Qc <- direct_sum(Qc, mat_one(n_ijkl[l], n_ijkl[l]) / n_ijkl[l])
          Time <- c(Time, rep(l, n_ijkl[l]))
        }
        n_ijkl <- n_ijkl[-(1 : l)]
        n_ijk. <- c(n_ijk., nijk.)
        Dose <- c(Dose, rep(k,nijk.))
        x2 <-  direct_sum(x2,  matrix(rep(1, nijk.), ncol = 1))
        qGama1 <- direct_sum(qGama1, matrix(rep(1,  nijk.), ncol = 1) / nijk.)
        Qb <-  direct_sum(Qb, mat_one(nijk., nijk.) / nijk.)
        qDelta2 <- direct_sum(qDelta2, mat_one(nijk.,  d_ijk[k]) / nijk.)
        nij.. <- nij.. + nijk.
      }
      n_ij.. <- c(n_ij.., nij..)
      x1 <- direct_sum(x1, matrix(rep(1, nij..), ncol = 1))
      qBeta1 <- direct_sum(qBeta1, matrix(rep(1, nij..), ncol = 1) / nij..)
      Qa <- direct_sum(Qa, mat_one(nij.., nij..) / nij..)
      qGama2 <- direct_sum(qGama2, mat_one(nij.., c_ij[j]) / nij..)
      Compound <- c(Compound, rep(j,nij..))
      d_ijk <- d_ijk[-(1 : k)]
      ni... <- ni... + nij..
    }
  QA <- Qa - mat_one(N, N) / N
  QB <- Qb - Qa
  QC <- Qc - Qb
  R <- diag(N) - Qc
  QBeta <- qBeta1 - mat_one(N, b_i) / N
  QGama <- qGama1 - qGama2
  QDelta <- qDelta1 - qDelta2
  retn <- list(
    X1 = x1,
    X2 = x2,
    X3 = x3,
    QBeta = QBeta,
    QGama = QGama,
    QDelta = QDelta
  )
  structure(retn, class = "NAMM")
}

# Function to get interaction matrix
get_coef <- function(chemicals, probes, expr_data, attr_data) {
  lev_str <- level_str2(chemicals, attr_data = attr_data)
  cc <- get_matrix(lev_str)
  ck_data <- cnr_data(
    chemicals,
    expr_data = expr_data,
    attr_data = attr_data,
    probes = probes,
    multicore = FALSE,
    store = FALSE)
  Y <- ck_data$expr_data
  X <- ck_data$attr_data
  # rownames(Y) <- get_genename(rownames(Y), get_organism(X))
  com_name <- unique(X$compound_abbr)
  doses <- substr(X$dose_level, start=1, stop=1)
  com_dose <- unique(glue("{X$compound_abbr} ({doses})"))
  timepoints <- gsub('([0-9]+) ([a-zA-Z])[a-zA-Z]*', '\\1\\2', X$time_level)
  com_dose_time <- unique(glue("{X$compound_abbr} ({doses}-{timepoints})"))
  coef_com <- Y %*% cc$QBeta
  colnames(coef_com) <- com_name
  coef_dose <- Y %*% cc$QGama
  colnames(coef_dose) <- com_dose
  coef_time <- Y %*% cc$QDelta
  colnames(coef_time) <- com_dose_time
  sam_avg <- avg_space(chemicals,
                       expr_data = Y,
                       attr_data = X,
                       probes = NULL,
                       dose = NULL,
                       time = NULL)
  return(list(compound = coef_com, dose = coef_dose, time = coef_time, samples = sam_avg$expr_data))
}

get_group <- function(comps_gr, attr_data) {

}






