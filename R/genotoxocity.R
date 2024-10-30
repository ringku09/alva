
# download_rawdata(comp_name="ET", data_type="R1LR", output_dir = "E:/NAAM/data_repeat/")
# download_rawdata(comp_name="CMA", data_type="R1LR", output_dir = "E:/NAAM/data_repeat/")

# files <- list.files("E:/NAAM/data_repeat/", full.names = TRUE)
#
# aaa <- process_compaunds(files,
#                          fc = TRUE,
#                          method = "rma",
#                          multicore = FALSE,
#                          store = TRUE,
#                          output_dir = "E:/NAAM/data_repeat")
#   dim(aaa$expression)
library(tidyverse)
library(rlang)
library(cli)
library(glue)
library(Matrix)
# install.packages("tidyverse")
# install.packages("RcppArmadillo")
# install.packages("RcppEigen")
# BiocManager::install("STRINGdb")

#options(timeout = 1000)
Rcpp::sourceCpp("mat_mult.cpp")
`%**%` <- function(x, y) matrixProd(x, y)

source("E:/lavana/NAAM/naam/R/1_down_data2.R")
source("E:/lavana/NAAM/naam/R/2_cnr_data.R")
source("E:/lavana/NAAM/naam/R/3_level_list.R")
source("E:/lavana/NAAM/naam/R/4_select_genes4.R")
source("E:/lavana/NAAM/naam/R/utils.R")
source("E:/lavana/NAAM/naam/R/5_probes2genes.R")
source("E:/lavana/NAAM/naam/R/get_naam_matrix.R")
source("E:/lavana/NAAM/naam/R/space_data.R")
source("E:/lavana/NAAM/naam/R/ggtheme.R")
source("E:/lavana/NAAM/naam/R/gene_enrich2.R")
source("E:/lavana/NAAM/naam/R/8_ppi3.R")
source("E:/lavana/NAAM/naam/R/heatmap-gene.R")
source("E:/lavana/NAAM/naam/R/biplot.R")
source("E:/lavana/NAAM/naam/R/plot_hottelingPCA.R")

  xx <- readr::read_csv('E:/lavana/NAAM/data_repeat/metadata.csv')
  yy <- readr::read_csv('E:/lavana/NAAM/data_repeat/expression.csv')
  expr_data <- df2matrix(dplyr::select(yy,-probes), dplyr::pull(yy, probes))
  attr_data <- xx
dim(expr_data)
ngeno <- c("AAA", "CBZ", "CCL4", "CFB", "EE", "FFB", "GFZ", "MCT", "MP", "MTS", "PB", "TAA", "WY")
nhepto <- c("3-MC", "AA", "AM", "ANIT", "APAP", "ASA", "IBU", "NIF", "PhB", "PPL", "RIF")
geno <- c("2NF", "AFB1", "AZP", "CBP","CPA","CSP")

# For single dose study
comps_gr <- list(ngeno =c("CBZ", "CCL4", "EE"),
                 nhepto = c("RIF", "AM", "IBU"),
                 geno = c("AZP", "CBP","CSP"))

# For repeat dose study
# comps_gr <- list(ngeno =c("TAA", "MP", "ET", "EE"),               #, "CCL4"
#                  nhepto = c("APAP", "AM", "ASA", "NIF"),      #, "NIF"
#                  geno = c("AZP", "CBP","CSP", "CPA"))          #, "CPA"

# comps_gr <- list(NGTH =c("TAA", "EE", "MP"),               #, "CCL4"
#                  GTH = c("NIF", "APAP", "AM"),      #, "NIF"
#                  NHC = c("AZP", "CBP","CPA"))          #, "CSP"


# comps_gr <- list(NGTH =c("TAA", "ET", "MP"),               #, "CCL4"
#                  GTH = c("ASA", "APAP", "AM"),      #, "NIF"
#                  NHC = c("AZP", "CBP","CPA"))          #, "CPA"

# Final set
comps_gr <- list(NGTH =c("TAA", "CCL4", "MP"),               #, "CCL4"
                 GTH = c("ASA", "APAP", "AM"),      #, "NIF"
                 NHC = c("AZP", "CBP","CPA"))          #, "CPA"


comps_gr <- list(ngeno= c("CCL4", "TAA", "ET","MP", "EE", "MCT"),
                 geno=c( "AZP", "CBP","CPA","CSP"))
system.time(
  res3 <- tgx_genes(comps_gr,
                    expr_data = expr_data,
                    attr_data = attr_data,
                    pval_cut = 0.05,
                    gr_diff = TRUE,
                    multicore = FALSE,
                    store = FALSE,
                    output_dir = missing_arg(),
                    error_call = caller_env())
)

res4 <- res3 %>%
  filter(sig_type == "DE")

print(res4, n=50)




