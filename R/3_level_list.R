
#' Structure of levels in compound data
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function `comp_str()` is used to find total number of levels in dose, time point within dose
#' and number of replications used in each dose-time combinations.
#'
#' `comp_str()` calculate total number of dose levels of a compound used for the experiment and
#' total number of time points for each doses, also total number of replication for each time
#' points within the dose. For instance, if a compound has 3 doses (Low, Medium, High) then dose
#' has 3 levels, if each dose level use 4 time points (3hr, 6hr, 9hr and 24hr, for single dose data)
#' then the time has levels (4, 4, 4) for 3 dose levels and finally if each dose-time combination
#' use 3 replications then replication has 3*4 or 4+4+4=12 levels (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3).
#' Therefore, the total of levels in each stage is a multiplication of earlier levels.
#'
#' @param comp Compound name or abbreviation
#' @param comp_attr Attribute data of the compound
#'
#' @return
#' A list of size for dose, time and replication
#' @export
#'
#' @examples
#' \donttest{
#' attr_data <- readr::read_csv("D:/project1_TGx/data_processing/data_dict.csv")
#' comp_str("APAP", attr_data)
#' }
comp_str <- function(comp, attr_data) {
  comp <- comp_abbr(comp)
  attr_data2 <- attr_data %>% filter(.data$compound_abbr == comp)
  dose_df <- attr_data2 %>%
    group_by(.data$dose_level) %>%
    nest()
  time_df <-lapply(dose_df$data, function(x) x %>% group_by(.data$time_level) %>% nest())
  dose_lev <- nrow(dose_df)
  time_lev <- sapply(time_df, nrow)
  rep_lev <- unlist(lapply(time_df, function(x) sapply(x$data, nrow)))
  lev_list <- list(dose_lev, time_lev, rep_lev)
  names(lev_list) <- c("dose", "time", "replication")
  return(lev_list)
}

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
expr_str <- function(..., attr_data, error_call = caller_env()) {
  comps_gr <- list2(...)
  if (length(comps_gr) == 1 && is_bare_list(comps_gr[[1]])) {
    comps_gr <- comps_gr[[1]]
  }
  comps_gr <- lapply(comps_gr, comp_abbr, error_call = error_call)
  group <- length(comps_gr)
  comp_gr <- unlist(lapply(comps_gr, length))
  comps_list <- as.list(unlist(comps_gr))
  zz <- lapply(comps_list, comp_str, attr_data = attr_data)
  dose_lev <- do.call(c,lapply(zz, "[[", 1))
  time_lev <- do.call(c,lapply(zz, "[[", 2))
  rep_lev <- do.call(c,lapply(zz, "[[", 3))
  comps_lev <- list(group, comp_gr, dose_lev, time_lev, rep_lev)
  names(comps_lev) <- c("group", "compound", "dose", "time", "replication")
  structure(comps_lev, class = "tggates")
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
design_matrix <- function(lev_str) {
  if (!is.list(lev_str) & !inherits(lev_str, "tggates")) {
    cli_abort(c("{.var lev_str} must be object of class `tggates` or `list`.",
                "x" = "The class {style_bold(col_cyan(backtick(class(lev_str))))} of \\
                {.var lev_str} is not supported.",
                "i" = "Please provide {.var lev_str} as a list of size 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."))
  } else if (length(lev_str)  !=  5) {
    cli_abort(c("The length of {.var lev_str} must be 5.",
                "i" = "You have supplied a list {.var lev_str} of size {length(lev_str)}, please
                make sure {.var lev_str} has a length of 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."))
  }
  a <- lev_str[[1]]
  b_i <- lev_str[[2]]
  c_ij <- lev_str[[3]]
  d_ijk <- lev_str[[4]]
  n_ijkl <- lev_str[[5]]
  b. <- sum(lev_str[[2]])
  c.. <- sum(lev_str[[3]])
  d... <- sum(lev_str[[4]])
  N <-  sum(lev_str[[5]])
  Group = Compound = Dose = Time = n_ijk. = n_ij.. = n_i... <- NULL
  x1 = x2 = x3 = x4 =  Qa = Qb = Qc = Qd <- matrix(nrow=0, ncol=0)
  qAlfa1 = qBeta1 = qBeta2 = qGama1 = qGama2 = qDelta1 = qDelta2 <-  matrix(nrow=0, ncol=0)
  ii <- 0
  for (i in 1 : a) {
    ni... <- 0
    for (j in 1 : b_i[i]) {
      nij.. <- 0
      for (k in 1 : c_ij[j]) {
        nijk. <- 0
        for (l in 1 : d_ijk[k]) {
          nijk. <- nijk. + n_ijkl[l]
          x4 <- direct_sum(x4, matrix(rep(1, n_ijkl[l]), ncol = 1))
          qDelta1 <- direct_sum(qDelta1, matrix(rep(1, n_ijkl[l]), ncol = 1) / n_ijkl[l])
          Qd <- direct_sum(Qd, mat_one(n_ijkl[l], n_ijkl[l]) / n_ijkl[l])
          Time <- c(Time, rep(l, n_ijkl[l]))
        }
        n_ijkl <- n_ijkl[-(1 : l)]
        n_ijk. <- c(n_ijk., nijk.)
        Dose <- c(Dose, rep(k,nijk.))
        x3 <-  direct_sum(x3,  matrix(rep(1, nijk.), ncol = 1))
        qGama1 <- direct_sum(qGama1, matrix(rep(1,  nijk.), ncol = 1) / nijk.)
        Qc <-  direct_sum(Qc, mat_one(nijk., nijk.) / nijk.)
        qDelta2 <- direct_sum(qDelta2, mat_one(nijk.,  d_ijk[k]) / nijk.)
        nij.. <- nij.. + nijk.
      }
      n_ij.. <- c(n_ij.., nij..)
      x2 <- direct_sum(x2, matrix(rep(1, nij..), ncol = 1))
      qBeta1 <- direct_sum(qBeta1, matrix(rep(1, nij..), ncol = 1) / nij..)
      Qb <- direct_sum(Qb, mat_one(nij.., nij..) / nij..)
      qGama2 <- direct_sum(qGama2, mat_one(nij.., c_ij[j]) / nij..)
      Compound <- c(Compound, rep(j+ii,nij..))
      d_ijk <- d_ijk[-(1 : k)]
      ni... <- ni... + nij..
    }
    n_i... <- c(n_i..., ni...)
    x1 <- direct_sum(x1, matrix(rep(1, ni...), ncol = 1))
    qAlfa1 <- direct_sum(qAlfa1, matrix(rep(1, ni...), ncol = 1) / ni...)
    Qa <- direct_sum(Qa, mat_one(ni..., ni...) / ni...)
    qBeta2 <- direct_sum(qBeta2, mat_one(ni..., b_i[i]) / ni...)
    c_ij  <- c_ij[-(1 : j)]
    Group <- c(Group, rep(i, ni...))
    ii <- ii + b_i[i]
  }
  QA <- Qa - mat_one(N, N) / N
  QB <- Qb - Qa
  QC <- Qc - Qb
  QD <- Qd - Qc
  R <- diag(N) - Qd
  QAlfa <- qAlfa1 - mat_one(N, a) / N
  QBeta <- qBeta1 - qBeta2
  QGama <- qGama1 - qGama2
  QDelta <- qDelta1 - qDelta2
  retn <- list(
    a = a,
    b. = b.,
    c.. = c..,
    d... = d...,
    N = N,
    n_i... = n_i...,
    n_ij.. = n_ij..,
    n_ijk. = n_ijk.,
    n_ijkl = n_ijkl,
    X1 = x1,
    X2 = x2,
    X3 = x3,
    X4 = x4,
    QA = QA,
    QB = QB,
    QC = QC,
    QD = QD,
    R = R,
    QAlfa = QAlfa,
    qAlfa = qAlfa1,
    QBeta = QBeta,
    qBeta = qBeta1,
    QGama = QGama,
    qGama = qGama1,
    QDelta = QDelta,
    qDelta = qDelta1,
    Group = Group,
    Compound = Compound,
    Dose = Dose,
    Time = Time
  )
  structure(retn, class = "tggates")
}
