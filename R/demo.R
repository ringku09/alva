library(tidyverse)
library(rlang)
library(cli)
library(glue)
# install.packages("tidyverse")
# install.packages("RcppArmadillo")
# install.packages("RcppEigen")
# BiocManager::install("STRINGdb")

#options(timeout = 1000)
Rcpp::sourceCpp("mat_mult.cpp")
`%**%` <- function(x, y) matrixProd(x, y)

source("E:/TGGATESlavana/R/1_down_data2.R")
source("E:/TGGATES/R/2_cnr_data.R")
source("E:/TGGATES/R/3_level_list.R")
source("D:/TGGATES/R/4_select_genes4.R")
source("E:/TGGATES/R/5_probes2genes.R")
source("E:/TGGATES/R/6_parallel.R")
source("E:/TGGATES/R/7_get_avgdata.R")
source("E:/TGGATES/R/utils.R")
source("E:/TGGATES/R/plot_genebar3.R")
source("E:/TGGATES/R/heatmap-gene.R")
source("E:/TGGATES/R/8_ppi3.R")
source("E:/TGGATES/R/gene_enrich2.R")
source("E:/TGGATES/R/volcano_plot.R")
source("E:/TGGATES/R/space_data.R")
source("E:/TGGATES/R/cordiagm.R")
source("E:/TGGATES/R/go_plot.R")
source("E:/TGGATES/R/get_sankey.R")
source("E:/TGGATES/R/ggtheme.R")
xx <- readr::read_csv('E:/rma_data/metadata.csv')
yy <- readr::read_csv('E:/rma_data/expression.csv')
#xx <- readr::read_csv('E:/GSH_realExpr/metadata.csv')
#yy <- readr::read_csv('E:/GSH_realExpr/expression.csv')
# D:/project1_TGx/rma_data/
expr_data <- df2matrix(dplyr::select(yy,-probes), dplyr::pull(yy, probes))
#expr_data <- expr_data[, -(1:108)]
attr_data <- xx

# comps_gr <- list(GHS = c("APAP", "MP", "NFZ", "BBZ", "CMA"),
#                  NOR = c("EME", "GBC", "HCB","INAH", "PH"))

comps_gr <- list(Positive = c("APAP", "MP", "NFZ", "BBZ", "CMA", "PHO","DEM"),
                 Negative = c("EME", "GBC", "HCB","INAH", "PH","GMC", "PEN"))
comps_gr <- list(GHSC = c("PHO", "DEM", "BSO", "BBZ"),
                 GHS = c("EME", "GBC", "HCB","INAH", "PH","GMC", "PEN"))
comps_gr <- list(GHS = c("APAP", "MP", "NFZ", "BBZ", "CMA"),
                 SYN = c( "BSO", "ET", "TAA"),
                 NOR = c("EME", "GBC", "HCB","INAH", "PH","GMC", "PEN"))

comps_gr <- list(GHSC = c("APAP", "MP", "NFZ"),
                 GHS = c("BSO", "ET", "TAA"),
                 NOR = c("GBC", "EME", "INAH"))

comps_gr <- list(Positive = c("APAP", "MP", "NFZ", "BBZ", "CMA", "PHO","DEM",
                              "MTZ", "AFB1", "TAA", "AZP", "PHE"),
                 Negative = c("EME", "GBC", "HCB","INAH", "PH","GMC", "PEN"))

a = 1:10
v = c("DDD", "gfg")
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
table(res3$sig_type)
res4 <- res3 %>%
    filter(sig_type == "DE")
res4$gene_name
  res4 <- res3 %>%
    filter(abs(`GHS-NOR`) >=1.3 | abs(`SYN-NOR`) >= 0.3 & sig_type == "DE")
mas5_gene <- res4$gene_name
rma_gene <- res4$gene_name

mas5_gene[!mas5_gene %in% rma_gene]
rma_gene[!rma_gene %in% mas5_gene]

filess <- list.files("E:/GSH_rawdata", full.names = TRUE)



process_compaunds(
  filess,
  fc = FALSE,
  method = "rma",
  multicore = TRUE,
  store = TRUE,
  output_dir = "E:/GSH_realExpr/"
)


# "fastgreedy", "walktrap", "edge.betweenness"
net_data <- get_netdata(
  res4,
  gene_col = "gene_name",
  p.value_col = "group",
  organism = "rat",
  score_threshold = 200,
  cluster_method = "edge.betweenness",
  version = "11.5")

pp <- plot_network(net_data, plot_layout = "dh", rm_noedge = TRUE)
ggsave(filename = "figures/network.png", plot = pp, width = 22, height = 16, dpi = 300, units = "cm")

table(net_data$vertices$gene_class)


hub_data <- get_hubdata(net_data, condition = "degree >= 30")
ppp <- plot_network(hub_data, plot_layout = "circle", rm_noedge = FALSE)
ggsave(filename = "figures/hubnet.png", plot = ppp, width = 11, height = 11, dpi = 300, units = "cm")

net_data$vertices$name[net_data$vertices$gene_class == "Cluster-2"]
xx <- cluster_enrich(res3, net_data, cluster = 4, path_n = 10)
plot_enrich(xx)
dev.off()
res4 %>%
  arrange(desc(`Positive-Negative`))
#  "Process", "Component", "Function", "KEGG", "Pfam", "InterPro" "WikiPathways" "
# All, Process, Component, Function, Keyword, KEGG, RCTM, Pfam, SMART, InterPro
# PMID
res5 <- res4 %>%
  filter(probe_id %in%  all_probe)
enrc_data <- get_enrich(gene_df = res4,
                        gene_col = "gene_name",
                        pvalue_col = "group",
                        category = "KEGG",
                        organism = "rat",
                        path_n = 20,
                        score_threshold = 200,
                        version = "11.5")

enr <- plot_enrich(enrc_data, plot_layout = "fr")

ggsave(filename = "figures/kegg_net2.png", plot = enr, width = 22, height = 14, dpi = 300, units = "cm")


enrich_cl2 <- cluster_enrich(gene_df = res4,
                           net_data = net_data,
                           cluster = 2,
                           gene_col = "gene_name",
                           pvalue_col = "group",
                           category = "Function",
                           organism = "rat",
                           path_n = 20,
                           score_threshold = 200,
                           version = "11.5")

#png("figures/enrich.png", width = 18, height = 12, res = 300, units = "cm")
#dev.off()

# data display boxplot
ck_data <- cnr_data(
  comps_gr,
  expr_data = expr_data,
  attr_data = attr_data,
  probes =all_probe,
  multicore = FALSE,
  output_dir = tempdir(),
  store = FALSE
)

ck_data <- get_space(comps_gr,
                     expr_data = expr_data,
                     attr_data = attr_data,
                     probes = res4$probe_id,
                     dose = "High",
                     time = "9 hr")

md <- ck_data$attr_data %>%
  mutate(sam_id = paste(compound_abbr,dose_level,time_level,sample_id, sep = "_"))
ge <- ck_data$expr_data
# colnames(ge) <- md$sam_id
xx <- calculate_silhouette_df(t(ge), md$com_class)
mean(xx$Silhouette_Width)
plot_silhouette_df(xx)
nrow(xx)

ge <- ge[,!colnames(ge) %in% xx$Sample_ID[which(xx$Silhouette_Width < -0.3)]]
md <- md %>%
  filter(!md$barcode %in% xx$Sample_ID[which(xx$Silhouette_Width < -0.3)])

print(md, n=350)

dim(GE)
pc <- prcomp(GE,
             center = TRUE,
             scale. = TRUE)
#library(factoextra)

fviz_eig(pc)
biplot(pc)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


meta <- ck_data$attribute
levstr <- expr_str(comps_gr, attr_data = meta)
desmat <- design_matrix(levstr)
length(ge["1368121_at",])
df <- data.frame(name = as.factor(desmat$Group), value = ge["1368121_at",],
                 com = as.factor(desmat$Compound))
library(viridis)
df %>%
  ggplot( aes(x=com, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(x=com, y=value, fill=name), size=1.5, alpha=0.5) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")+
  facet_wrap(~name)

ge2 <- get_simdata(n_gene=1,n_com = c(5,5))[[1]]
df <- data.frame(name = as.factor(rep(c("A","B"),each=180)), value = as.vector(ge2))
#================================ Test compounds  =============================#
 pos <- c("CSP", "DOX", "DFNa", "NPX", "CPX", "AFB1","AA",
          "AM", "CPA", "CCL4", "MTZ", "AZP")
 neg <- c("ACZ","APL","CBZ", "DEX", "DZP", "FAM", "IM",
          "NIF", "RAN", "AMT", "CAP", "CIM","CFB", "PHE")
 comps_gr <- list(Positive = c("APAP", "MP", "NFZ", "BBZ", "CMA", "PHO","DEM",
                               "ETN", "VPA", "ASA", "MFM", "OPZ", "FLX",
                               "CSP", "DOX", "DFNa", "NPX", "CPX", "AFB1","AA",
                               "CPA", "AM", "CCL4", "MTZ", "AZP", "TAA"),
                  Negative = c("EME", "GBC", "HCB","INAH", "PH","GMC", "PEN",
                               "ACZ","APL","CBZ", "DEX", "DZP", "FAM", "IM",
                               "NIF", "RAN", "AMT", "CAP", "CIM","CFB", "PHE"))

 comps_gr <- list(Positive = c("APAP", "MP", "NFZ", "BBZ", "CMA", "PHO","DEM",
                               "MTZ", "PHE", "NIF", "CBZ", "AFB1", "TAA", "OPZ"),
                  Negative = c("EME", "GBC", "HCB","INAH", "PH","GMC", "PEN",
                               "APL", "ACZ", "FAM", "CAP", "DEX", "CFB", "DZP", "CIM"))

trts_data <- expand_data(unlist(comps_gr),
                        expr_data=expr_data,
                        attr_data=attr_data,
                        multicore = TRUE,
                        store = TRUE,
                        output_dir = "D:/project1_TGx/rma_test",
                        error_call = caller_env())

files <- list.files("D:/project1_TGx/rma_test/", full.names = TRUE)
tgx_data <- process_compaunds(files = files[17],
                              store = TRUE,
                              output_dir = "D:/project1_TGx/rma_test")
dim(tgx_data$expression)
dim(tgx_data$metadata)
unique(attr_data$compound_name)
expn2 <- cbind(expr_data, tgx_data$expression)
meta2 <- rbind(attr_data, tgx_data$metadata)

expn2 %>%
  tibble::as_tibble(rownames = "probes") %>%
  readr::write_csv(file = paste(output_dir, "expression2.csv", sep = "/"))
meta2 %>%
  readr::write_csv(file = paste(output_dir, "metadata2.csv", sep = "/"))
#======================================  package development ==========================


library(clusterProfiler)
hsa_kegg = download_KEGG("rno")









usethis::use_r("simulation.R")
usethis::use_package("doParallel")  #, type = "Suggests"
devtools::document()
devtools::build()
devtools::check()


#compds %in% compounds()$COMPOUND_ABBREVIATION
compds <-c("2NF", "NMOR")
out <- 'D:/project1_TGx/data_processing'
# dt <- lapply(compds, download_rawdata, data_type = "R1LS", output_dir = out)
#list.files(out, full.names = TRUE)
tgx_data <- get_TGGATES(compds, "rat", "R1LS", log2FC = TRUE,.parallel=3, out)
head(tgx_data$ge_data)
head(tgx_data$data_dict)
str(tgx_data$data_dict)



dd <- cnr_data(comps_gr, expr_data = expr_data, attr_data = attr_data)
dd$expression

newdt <- cnr_data(
  comps_gr,
  expr_data = getdt$expression,
  attr_data = getdt$attribute,
  multicore = 3,
  output_dir = "D:/project1_TGx/data_processing5/",
  store = TRUE,
  error_call = caller_env()
)

check_data <- cnr_data(
  comps_gr,
  expr_data = newdt$expression,
  attr_data = newdt$attribute,
  multicore = 3,
  output_dir = tempdir(),
  store = FALSE
)

unique(newdt$attribute$compound_abbr)


#install.packages("devtools")
#devtools::install_github("hadley/devtools")

library(devtools)
packageVersion("devtools")

file.exists("~/.ssh/id_rsa.pub")
usethis::use_mit_license()
usethis::create_package("C:/Users/RINGKU/Desktop/TGGATES")
usethis::use_github()
gh_toke=n_help()
usethis::use_readme_rmd()
usethis::use_r("reexport.R")
file.exists("C:/Users/RINGKU/.ssh/id_rsa.pub")
usethis::edit_r_environ()
usethis::use_package("zip")


###
set.seed(100)
x1 = rnorm(100000)
x2 = rnorm(100000)
x3 = rnorm(100000)
x4 = rnorm(100000)
x5 = rnorm(100000)
x6 = rnorm(100000)
x7 = rnorm(100000)
x8 = rnorm(100000)
x9 = rnorm(100000)
x10 = rnorm(100000)
y = 2*x1+0.5*x2+x3+1.5*x5+x7*x8+2*x9+0.05*x10
raw_data <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, y)
#head(raw_data)
#summary(lm(y~., data = raw_data))
meta_data <- read.csv("D:/project1_TGx/data_processing/2-nitrofluorene.Rat.in_vivo.Liver.Single/Attribute.tsv")
head(meta_data)


usethis::use_data(comp_data, compress = "xz")
R.version.string

usethis::create_github_token()
gitcreds::gitcreds_set()

usethis::use_vignette("my-vignette")

library(tidysum)

data("raw_data")
head(raw_data)

#getOption('timeout')
#options(timeout=1000)
fileName <- "ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Single/N-nitrosomorpholine.Rat.in_vivo.Liver.Single.zip"
temp <- "D:\\temp_download\\2-nitrofluorene.Rat.in_vivo.Liver.Single.zip"
download.file(fileName,temp)
unzip("D:\temp_download\expr_data", exdir = "D:\\temp_download")
file.remove(list.files(td, full.names = TRUE))
raw_path[1]
files(data_raw)
unlink(temp)
library(plyr)
install.packages("plyr")

usethis::use_pipe()
usethis::use_package("dplyr")
devtools::document()
devtools::build()
devtools::check()

plot(res4$`GHS_C-GHS`, res4$`GHS-Non_GHS`)
library(TGGATES)
compounds()
compound_link("2NF", "R1LS")
?download_rawdata()
temp_link <- compound_link("2NF", "R1LS")
rw_dt <- download_rawdata(temp_link, "D:/project1_TGx/data_processing")

chem_data <- read.table(file = "D:/project1_TGx/data_processing/2-nitrofluorene.Rat.in_vivo.Liver.Single/Attribute.tsv", sep = '\t', header = TRUE)
names(chem_data)
head(chem_data)
chem_data[,c("BARCODE","COMPOUND_NAME", "COMPOUND.Abbr.","DOSE_LEVEL","SACRI_PERIOD","INDIVIDUAL_ID")]
temp_link <- compound_link("2NF", "R1LS")
cel_path <- download_rawdata(temp_link, getwd())
cel_expr <- cel_to_exprs(cel_path$CEL_files[1:2], data_type = "R1LS")
head(cel_expr)
?cel_to_exprs
