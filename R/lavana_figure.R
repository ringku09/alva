
library(cowplot)
library(factoextra)
library(ggVennDiagram)
library(ComplexUpset)
library(pheatmap)
#---------------------------  Figure 1  -------------------------------------------
mod_alfa05 <- tgx_genes(comps_gr,
                  expr_data = expr_data,
                  attr_data = attr_data,
                  pval_cut = 0.05)
mod_alfa1e5 <- tgx_genes(comps_gr,
                  expr_data = expr_data,
                  attr_data = attr_data,
                  pval_cut = 1e-5)
gene_05 <- mod_alfa05 %>%
  filter(sig_type == "DE")
gene_1e5 <- mod_alfa1e5 %>%
  filter(sig_type == "DE")
naov_data1 <- get_coef(unlist(comps_gr),
               probes = mod_alfa05$probe_id,
               expr_data = expr_data,
               attr_data = attr_data)
naov_data2 <- get_coef(unlist(comps_gr),
                      probes = gene_05$probe_id,
                      expr_data = expr_data,
                      attr_data = attr_data)
naov_data3 <- get_coef(unlist(comps_gr),
                      probes = gene_1e5$probe_id,
                      expr_data = expr_data,
                      attr_data = attr_data)

################  n component
df <- data.frame(naov_data3$compound)
nBartlett(df, nrow(df))
rndLambdaF(df)


spca <- SamplePCA(t(df))
lambda <- spca@variances[1:(ncol(df)-1)]
bsDimension(lambda)

ag.obj <- AuerGervini(spca)
agDimension(ag.obj)

ag.obj <- AuerGervini(spca)
agDimension(ag.obj)

f <- makeAgCpmFun("Exponential")
agfuns <- list(twice=agDimTwiceMean, specc=agDimSpectral,
                 km=agDimKmeans, km3=agDimKmeans3,
                 tt=agDimTtest, tt2=agDimTtest2,
                 cpt=agDimCPT, cpm=f)
compareAgDimMethods(ag.obj, agfuns) # compare the list of all criteria
pts <- screeplot(spca)
lines(pts, bs, type='b', col='blue', lwd=2, pch=16)
lines(pts[-NC], bs0, type='b', col='red', lwd=2, pch=16)

plot(ag.obj, agfuns)

######################################

pca1 <- prcomp(naov_data1$compound, center = TRUE, scale. = TRUE)
pca2 <- prcomp(naov_data2$compound, center = TRUE, scale. = TRUE)
pca3 <- prcomp(naov_data3$compound, center = TRUE, scale. = TRUE)
scree1 <- fviz_screeplot(pca1, addlabels = TRUE, barfill = "gray", barcolor ="black", linecolor ="#FF00FF",
                         n = which(cumsum(pca1$sdev^2/sum(pca1$sdev^2))>=0.9)[1],
                         ylim = c(0,65), title = "", ylab = "", xlab ="")
scree2 <- fviz_screeplot(pca2, addlabels = TRUE, barfill = "gray", barcolor ="black", linecolor ="#FF00FF",
                         n = which(cumsum(pca2$sdev^2/sum(pca2$sdev^2))>=0.9)[1],
                         ylim = c(0,65), title = "", ylab = "", xlab ="")
scree3 <- fviz_screeplot(pca3, addlabels = TRUE, barfill = "gray", barcolor ="black", linecolor ="#FF00FF",
                         n = which(cumsum(pca3$sdev^2/sum(pca3$sdev^2))>=0.9)[1],
                         ylim = c(0,65), title = "", ylab = "", xlab ="")


pca <- prcomp(naov_data1$compound, center = TRUE, scale. = TRUE)
com_expr <- unlist(lapply(str_split(colnames(naov_data1$compound), " "),"[",1))
com_group <- names(comps_gr)[get_class(comps_gr, com_expr)]

biplot1 <- bi_plot(comps_gr,naov_data1$compound)
biplot2 <- bi_plot(comps_gr,naov_data2$compound)
biplot3 <- bi_plot(comps_gr,naov_data3$compound)

fig1_a2c <- plot_grid(scree1,
                      scree2,
                      scree3,
                      ncol = 3, nrow = 1, labels = letters[1:3],  label_size = 20)
fig1_a2c <- plot_grid(arrangeGrob(fig1_a2c,
                               left = textGrob("Explained variances (%)", rot = 90, vjust = 1,
                                               gp = gpar(fontsize = 14, fontface = "bold")),
                               bottom = textGrob("Dimensions",
                                                 gp = gpar(fontsize = 14, fontface = "bold"))))

fig1_d2f <- plot_grid(biplot1 + theme(legend.position="none"),
                      biplot2 + theme(legend.position="none"),
                      biplot3 + theme(legend.position="none"),
                      ncol = 3, nrow = 1, labels = letters[4:6],  label_size = 20)
fig1_legend <- common_legend(biplot1)
fig1_d2f <- plot_grid(fig1_d2f,
                      fig1_legend, nrow = 2,rel_heights = c(0.92, 0.08))
fig1 <- plot_grid(fig1_a2c,fig1_d2f,nrow = 2)

ggsave(filename = "E:/lavana/results/Figure1.jpeg", plot = fig1, width = 190*2, height = 140*2, dpi = 300, units = "mm")



biplot_a <- bi_plot2(comps_gr,naov_data2$samples, dose = "L")
biplot_b <- bi_plot2(comps_gr,naov_data2$time, dose = "L")
biplot_c <- bi_plot2(comps_gr,naov_data2$dose, dose = "L")

biplot_d <- bi_plot2(comps_gr,naov_data2$samples, dose = "M")
biplot_e <- bi_plot2(comps_gr,naov_data2$time, dose = "M")
biplot_f <- bi_plot2(comps_gr,naov_data2$dose, dose = "M")

biplot_g <- bi_plot2(comps_gr,naov_data2$samples, dose = "H")
biplot_h <- bi_plot2(comps_gr,naov_data2$time, dose = "H")
biplot_i <- bi_plot2(comps_gr,naov_data2$dose, dose = "H")


fig2_a2f <- plot_grid(biplot_a+ theme(legend.position="none"),
                      biplot_b+ theme(legend.position="none"),
                      biplot_c+ theme(legend.position="none"),
                      biplot_d+ theme(legend.position="none"),
                      biplot_e+ theme(legend.position="none"),
                      biplot_f+ theme(legend.position="none"),
                      biplot_g+ theme(legend.position="none"),
                      biplot_h+ theme(legend.position="none"),
                      biplot_i+ theme(legend.position="none"),
                      ncol = 3, nrow = 3, labels = letters[1:9],  label_size = 20)
fig2_legend <- common_legend(biplot_a)
fig2 <- plot_grid(fig2_a2f,
                      fig2_legend, nrow = 2,rel_heights = c(0.98, 0.02))

ggsave(filename = "E:/lavana/results/Figure2.jpeg", plot = fig2, width = 190*2, height = 190*2, dpi = 300, units = "mm")



gene_05 <- mod_alfa05 %>%
  filter(sig_type == "DE")

naov_data2 <- get_coef(unlist(comps_gr),
                       probes = NULL,
                       expr_data = expr_data,
                       attr_data = attr_data)
fig3 <- plot_pcaellips(naov_data2$compound,res3, alpha = 0.00001, top=50, pts = 200)
ggsave(filename = "E:/lavana/results/Figure3.jpeg", plot = fig3, width = 190, height = 190, dpi = 300, units = "mm")



fig4.1 <- plot_pcamhat( pca_data = naov_data2$compound, gene_tab = gene_05, k=2, alpha = 0.05)
fig4.2 <- plot_pcamhat( pca_data = naov_data2$compound, gene_tab = gene_05, k=3, alpha = 0.05)
fig4.3 <- plot_pcamhat( pca_data = naov_data2$compound, gene_tab = gene_05, k=4, alpha = 0.05)
fig4.4 <- plot_pcamhat( pca_data = naov_data2$compound, gene_tab = gene_05, k=5, alpha = 0.05)


fig4_a2d <- plot_grid(fig4.1 + labs(x = "", y = ""),
                      fig4.2 + labs(x = "", y = ""),
                      fig4.3 + labs(x = "", y = ""),
                      fig4.4 + labs(x = "", y = ""),
                      ncol = 2, nrow = 2, labels = letters[1:4],  label_size = 20)
fig4_a2d <- plot_grid(arrangeGrob(fig4_a2d,
                                  left = textGrob(expression("Hotelling’s T"^2), rot = 90, vjust = 1,
                                                  gp = gpar(fontsize = 18, fontface = "bold")),
                                  bottom = textGrob(expression("Gene (ordered by -log"[10]*"P-value)"),
                                                    gp = gpar(fontsize = 18, fontface = "bold"))))

ggsave(filename = "E:/lavana/results/Figure4.jpeg", plot = fig4_a2d, width = 200*2, height = 180*2, dpi = 300, units = "mm")



pca2 <- get_pcaT2(pca_data = naov_data2$compound, gene_tab = gene_05, k=2, alpha = 0.05)
pca3 <- get_pcaT2(pca_data = naov_data2$compound, gene_tab = gene_05, k=3, alpha = 0.05)
pca4 <- get_pcaT2(pca_data = naov_data2$compound, gene_tab = gene_05, k=4, alpha = 0.05)
pca5 <- get_pcaT2(pca_data = naov_data2$compound, gene_tab = gene_05, k=5, alpha = 0.05)


gene2 <- pca2$pca_data$probe_id[pca2$pca_data$gene == "DE"]
gene3 <- pca3$pca_data$probe_id[pca3$pca_data$gene == "DE"]
gene4 <- pca4$pca_data$probe_id[pca4$pca_data$gene == "DE"]
gene5 <- pca5$pca_data$probe_id[pca5$pca_data$gene == "DE"]

length(gene2)

probe_set <- list(`K = 2` = pca2$pca_data$gene_symbol[pca2$pca_data$gene == "DE"],
                  `K = 3` = pca2$pca_data$gene_symbol[pca2$pca_data$gene == "DE"],
                  `K = 4` = pca2$pca_data$gene_symbol[pca2$pca_data$gene == "DE"],
                  `K = 5` = pca2$pca_data$gene_symbol[pca2$pca_data$gene == "DE"])
sankey_kgbp <- enrich_sankey(probe_set,
                              left_path = "KEGG",
                              right_path = "Process",
                              header = c("KEGG", "GO : BP"),
                              height = 900,
                              width = 2000,
                              font_size = 26,
                              header_size = 26)
htmlwidgets::saveWidget(sankey_kgbp, "E:/lavana/results/sankey.html")
webshot::webshot("E:/lavana/results/sankeyNetwork.html","E:/lavana/results/sankey.png", vwidth = 2000, vheight = 900)

#


all_genes <- Reduce(union, list(gene2, gene3, gene4, gene5))
upset_data <- tibble(
  genes = all_genes,
  `k = 2` = all_genes %in% gene2,
  `k = 3` = all_genes %in% gene3,
  `k = 4` = all_genes %in% gene4,
  `k = 5` = all_genes %in% gene5
)

upset_plot <- upset(
  upset_data,
  c("k = 2","k = 3", "k = 4", "k = 5"),
  name = "Number of PCs",
  base_annotations = list(
    'Intersection Size'=intersection_size(
      counts=TRUE
    )),
  set_sizes=(
    upset_set_size(
      geom = geom_bar(
        aes( x=group),
        width=0.5),position='right')))

veen_data <- list(
  `k = 2` = gene2,
  `k = 3` = gene3,
  `k = 4` = gene4,
  `k = 5` = gene5
)

venn_plot <- ggVennDiagram(veen_data, color = "black", lwd = 0.3, lty = 1, size= 0.5,
                           label_percent_digit = 1, label_size = rel(3), set_size = rel(4)) +
  scale_x_continuous(expand = expansion(mult = 0.12))+
  scale_fill_gradient(low = "#F4FAFE", high = "#FF4500")

fig5 <- ggplotify::as.grob(upset_plot + annotation_custom(ggplotGrob(venn_plot),
                                                            xmin = -5, ymin= -400,
                                                            xmax = 25))

ggsave(filename = "E:/lavana/results/Figure5.jpeg", plot = fig5, width = 190, height = 150, dpi = 300, units = "mm")

write.csv(pca2$pca_data %>% filter(gene == "DE"), "E:/lavana/results/pca_df.csv")


dd <- read.csv(file = "E:/lavana/results/gene_df.csv",header = TRUE)


df_ppi <- pca2$pca_data %>% filter(gene == "DE")
net_edgebet <- get_netdata(
  gene_df = df_ppi,
  probe_col = "probe_id",
  gene_col = "gene_symbol",
  p.value_col = "p_value",
  organism = "rat",
  score_threshold = 200,
  cluster_method = NULL,
  version = "12")

ppi_class <- dd$class[match(net_edgebet$vertices$probes, dd$probe_id)]
community <- as.numeric(gsub("[^0-9]", "", ppi_class))
net_edgebet$vertices$community <- community
net_edgebet$vertices$gene_class <- ppi_class

fig.5a <- plot_network(net_data = net_edgebet, plot_layout = "dh", rm_noedge = TRUE)

deg <- get_hubdata(net_edgebet, condition = "degree >= 26")
eig <- get_hubdata(net_edgebet, condition = "eigenes >= 0.5")
deg_plot <- plot_network(deg, plot_layout = "circle",arc = 0,
                       rm_noedge = FALSE, show_legend = FALSE)
eig_plot <- plot_network(eig, plot_layout = "sphere",arc = 0,
                       rm_noedge = FALSE, show_legend = FALSE)

fig.5bc <- plot_grid(deg_plot, eig_plot, labels = c('b', 'c'), label_size = 20, ncol = 1)
fig.5abc <- plot_grid(fig.5a, fig.5bc, labels = c('a', ''), rel_widths = c(0.7,0.3), label_size = 20, ncol = 2)

ggsave(filename = "E:/lavana/results/Figure12d.jpeg", plot = fig.5abc, width = 180*2, height = 200, dpi = 300, units = "mm")

write.csv(net_edgebet$vertices, "E:/lavana/results/network_df.csv")
#--------------------------   heatmap

fig_s1.7a <- tgx_heatmap(comps_gr,
                         expr_data = expr_data,
                         attr_data = attr_data,
                         probes = gene2,
                         dose =  "High",
                         time = NULL,
                         annotation_compound = TRUE,
                         annotation_legend = TRUE)


head(naov_data2$time)


#----------------------------------  plsa ------------------------------------------
remove_na <- function(x) {
  x[!is.na(x)]
}
naov_data2 <- get_coef(unlist(comps_gr),
                       probes = NULL,
                       expr_data = expr_data,
                       attr_data = attr_data)
gene2 <- fig3$data$probe_id[fig3$data$gene_symbol %in% remove_na(fig3$data$gene_symbol[fig3$data$gene == "DE"])]
vv <- naov_data2$compound[gene2, ]
rownames(vv) <- get_genename(rownames(vv), organism = "rat")
# GCmat <- round(exp(vv)*50)
GCmat <- round(1/(1+exp(-vv))*100) #vv_mat/sum(vv_mat)
 pp <- PHVM(GCmat, k = 2)

  # grey(seq(1, 0.3, length = 1000))
#-----------------------------------  cocluster plot -----------------------------
cc <- BioMCocls(GCmat = GCmat,
                CoCls = 2,
                GCjointProb = pp$GCjointProb,
                GclsMem = pp$Gclust,
                CclsMem = pp$Cclust)
pr_GC <- scales::rescale(cc$Co_cls_JointProb, to = c(0, 1))

pheatmap::pheatmap(t(pr_GC), cluster_row = FALSE, cluster_col = FALSE, border_color= NA,
                   color = grey(seq(1, 0, length = 1000)),
                   width = 15, height = 10, filename = "E:/lavana/results/Figure6.jpeg")
#---------------------------  gene and compound cluster ------------------------------
gene_member <- sapply(strsplit(rownames(as.matrix(unlist(pp$Gclust))),".",fixed=TRUE),function(x)x[2])
gg_class <- pp$GzProb[gene_member,]
pheatmap::pheatmap(gg_class, cluster_row = FALSE, cluster_col = FALSE, border_color= NA,
                   color = grey(seq(1, 0, length = 1000)),
                   width = 10, height = 10, filename = "E:/lavana/results/Figure9.jpeg")

com_member <- sapply(strsplit(rownames(as.matrix(unlist(pp$Cclust))),".",fixed=TRUE),function(x)x[2])
cc_class <- pp$CzProb[com_member,]
pheatmap::pheatmap(cc_class, cluster_row = FALSE, cluster_col = FALSE, border_color= NA,
                   color = grey(seq(1, 0, length = 1000)),
                   width = 10, height = 10, filename = "E:/lavana/results/Figure10.jpeg")

#------------------------------------  top table  --------------------------------------------------
long_df <- as.data.frame(pr_GC) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Compound", values_to = "Score") %>%
  arrange(desc(Score))
top_values <- long_df %>%
  group_by(Compound) %>%
  arrange(desc(Score)) %>%
  slice_head(n = 10) %>%
  ungroup()  # Remove grouping
# write.csv(top_values, "E:/lavana/results/top10.csv")
top_10_max_values <- long_df %>%
  group_by(Compound) %>%
  arrange(desc(Score)) %>%
  slice_head(n = 10) %>%
  ungroup()  # Remove grouping

print(df, n = 10)

df <- top_10_max_values %>%
  dplyr::group_by(Compound, Gene) %>%
  dplyr::summarise(value = n()) %>%
  dplyr::rename(from = Compound, to = Gene)
cordm <- net_enrichment(comps_gr, df,lab_size = 0.8)

ggsave(filename = "E:/lavana/results/Figure11.jpeg", plot = cordm, width = 180, height = 180, dpi = 300, units = "mm")

##------------------  Gene Table ------------------------------------
gene_maxprob <- apply(scales::rescale(pp$GzProb, to = c(0, 1)), 1, max)
gene_list <- lapply(pp$Gclust, function(x) names(x))
df_gene <- data.frame(gene_name = unlist(gene_list, use.names = FALSE))
df_gene$class <- rep(names(pp$Gclust), times = sapply(pp$Gclust, length))
df_gene$class_prob <- gene_maxprob[match(names(gene_maxprob), df_gene$gene_name)]

tox_df <- gene_05 %>%
  filter(probe_id %in% gene2) %>%
  select(c("probe_id", "entrez_id", "gene_name", "gene_fnane", "log10p"))

gene_df <- left_join(tox_df, df_gene, by = 'gene_name') %>%
  relocate(class_prob, .before = log10p)
write.csv(gene_df, "E:/lavana/results/gene_df.csv")

#-------------------------  Compound Table  -----------------------------
com_maxprob <- apply(scales::rescale(pp$CzProb, to = c(0, 1)), 1, max)

com_list <- lapply(pp$Cclust, function(x) names(x))
df_com <- data.frame(com_name = unlist(com_list, use.names = FALSE))
df_com$class <- rep(names(com_list), times = sapply(com_list, length))
df_com$class_prob <- com_maxprob[match(names(com_maxprob), df_com$com_name)]
write.csv(df_com, "E:/lavana/results/com_df.csv")
#----------------------     Probability data  -------------------------
write.csv(cc$Co_cls_JointProb, "E:/lavana/results/co_clust.csv")
write.csv(pp$GCjointProb, "E:/lavana/results/joint_prob.csv")
write.csv(pp$GzProb, "E:/lavana/results/gene_prob.csv")
write.csv(pp$CzProb, "E:/lavana/results/com_prob.csv")
#-----------------------------------------------------------



#----------------------------  pathway analysis  ----------------------------
print(fig2)
gene_df2 <- gene_05 %>%
  filter(probe_id %in% gene2)
dim(gene_df2)
enrc_kegg <- get_enrich(gene_df = gene_df2,
                        gene_col = "gene_name",
                        pvalue_col = "p_values",
                        category = "KEGG",
                        organism = "rat",
                        path_n = 15,
                        score_threshold = 200,
                        version = "12")

print(enrc_kegg$vertices, n=15)

fig.4a <- plot_enrich(enrc_kegg, plot_layout = "fr")

kegg_ovall <- enrich_genes(gene_df = gene_df2,
                         gene_col = "gene_name",
                         category = "KEGG",
                         organism = "rat",
                         path_n = 50,
                         score_threshold = 200,
                         version = "12")

process_ovall <- enrich_genes(gene_df = gene_df2,
                           gene_col = "gene_name",
                           category = "Process",
                           organism = "rat",
                           path_n = 100,
                           score_threshold = 200,
                           version = "12")

kegg_ovall$path_name <- split2_merge(kegg_ovall$path_name)
kegg <- ggplot(kegg_ovall[1:10,], aes(-log10(fdr), reorder(path_name, n_genes), fill = n_genes)) +
  geom_bar(stat = "identity")+
  labs(x=expression("-log"[10]*" (FDR):"), y= "KEGG",
       fill = "#Genes")+
  scale_fill_gradient(low = "royalblue", high = "royalblue4")+
  theme_Publication(tag_size = rel(1.5),
                    leg_size = rel(1),
                    #x_lab = element_blank(),
                    #y_lab = element_blank()
  )+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")

process_ovall$path_name <- split2_merge(process_ovall$path_name)
process <- ggplot(process_ovall[1:10,], aes(-log10(fdr), reorder(path_name, n_genes), fill = n_genes)) +
  geom_bar(stat = "identity")+
  labs(x=expression("-log"[10]*" (FDR):"), y= "GO : BP",
       fill = "#Genes")+
  scale_fill_gradient(low = "tomato", high = "tomato4")+
  theme_Publication(tag_size = rel(1.5),
                    leg_size = rel(1),
                    #x_lab = element_blank(),
                    #y_lab = element_blank()
  )+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")

kegg_proc <- plot_grid(kegg, process,ncol = 2, labels = c("a", "b"), label_size = 20)
ggsave(filename = "E:/lavana/results/Figure7.jpeg", plot = kegg_proc, width = 190*2, height = 220, dpi = 300, units = "mm")



wiki_ovall <- enrich_genes(gene_df = gene_df2,
                              gene_col = "gene_name",
                              category = "WikiPathways",
                              organism = "rat",
                              path_n = 100,
                              score_threshold = 200,
                              version = "12")
rctm_ovall <- enrich_genes(gene_df = gene_df2,
                           gene_col = "gene_name",
                           category = "RCTM",
                           organism = "rat",
                           path_n = 100,
                           score_threshold = 200,
                           version = "12")

wiki_ovall$path_name <- split2_merge(wiki_ovall$path_name)
wiki <- ggplot(wiki_ovall[1:6,], aes(-log10(fdr), reorder(path_name, n_genes), fill = n_genes)) +
  geom_bar(stat = "identity")+
  labs(x=expression("-log"[10]*" (FDR):"), y= "Wikipathways",
       fill = "#Genes")+
  scale_fill_gradient(low = "royalblue", high = "royalblue4")+
  theme_Publication(tag_size = rel(1.5),
                    leg_size = rel(1),
                    #x_lab = element_blank(),
                    #y_lab = element_blank()
  )+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")

rctm_ovall$path_name <- split2_merge(rctm_ovall$path_name)
rctm <- ggplot(rctm_ovall[1:10,], aes(-log10(fdr), reorder(path_name, n_genes), fill = n_genes)) +
  geom_bar(stat = "identity")+
  labs(x=expression("-log"[10]*" (FDR):"), y= "Reactome",
       fill = "#Genes")+
  scale_fill_gradient(low = "tomato", high = "tomato4")+
  theme_Publication(tag_size = rel(1.5),
                    leg_size = rel(1),
                    #x_lab = element_blank(),
                    #y_lab = element_blank()
  )+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")

rctm_wiki <- plot_grid(rctm, wiki,ncol = 2, labels = c("a", "b"), label_size = 20)
ggsave(filename = "E:/lavana/results/Figure8.jpeg", plot = rctm_wiki, width = 190*2, height = 220, dpi = 300, units = "mm")

write.csv(func_ovall, "E:/lavana/results/func_ovall.csv")
##---------------------------------------------------------------




