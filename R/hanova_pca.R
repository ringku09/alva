
library(factoextra)

cc <- cnr_data(comps_gr, probes = res4$probe_id,
               expr_data = expr_data, attr_data = attr_data)

bb <- get_coef(unlist(comps_gr),
               probes = res3$probe_id,
               expr_data = expr_data,
               attr_data = attr_data)

pca <- prcomp(bb$compound, center = TRUE, scale. = TRUE)
fviz_pca_var(pca, col.var = "black", axes = c(1,2))




fviz_screeplot(pca, addlabels = TRUE, n = 5, ylim = c(0, 60))

fviz_pca_biplot(pca, repel = FALSE, axes = c(1,2),
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

print(ck_data$attr_data,n=357)

yy <- t(ck_data$expression[res4$probe_id, ] )
xx <- cbind(design_mat$X1,design_mat$X2,design_mat$X3,design_mat$X4)
mod <- lm(yy~xx)
b <- coefficients(mod)
vv <- t(b[2:18,])

vv <- ck_data$expr_data %*% design_mat$QBeta
vv <- ck_data$expression[res4$probe_id, ] %*% design_mat$qBeta
dim(vv)
# vv <- vv[-c(6,9), ]
com <- unique(ck_data$attr_data$compound_abbr)
com_dose <- paste(unique(ck_data$attribute$compound_abbr),rep(unique(ck_data$attribute$dose_level), times =14))
colnames(vv) <- com
rownames(vv) <- res4$gene_name
# colnames(vv) <- paste(rep(unlist(comps_gr),each=3),rep(c("L","M","H"), times=12), sep = "_")
dim(vv)
cc <- net_data$vertices %>%
  filter(gene_class != "others-38")
net_data$vertices$probes
pca <- prcomp(vv, center = TRUE, scale. = TRUE)
fviz_eig(pca)
str(pca)
biplot(pca, scale = TRUE)

jbiplot(1, 2, vv,
        gnames = res4$gene_name,
        loadingnames = com)

fviz_pca_biplot(pca, repel = FALSE, axes = c(1,2),
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
dim(vv)
hist(vv[30, ])
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)


length(unique(colnames(vv)))
fviz_pca_biplot(pca,
                #col.ind = as.factor(vv2$cluster),
                palette = "jco",
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species")

fviz_pca_biplot(pca,
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                #fill.ind =as.factor(vv2$cluster),
                col.ind = "black",
                # Color variable by groups
                #col.var = rep(c("red","green","blue"),times=c(7,3,7)),

                legend.title = list(fill = "Species", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors


vv2 <- kmeans(x=vv, centers=3)
vv2$cluster[vv2$cluster == 2]
library(mclust)
mclustBIC(vv)
zz <- Mclust(vv, model = "VEE")
zz$classification
zz$classification[zz$classification == 1]


tt <- tight.clust(vv, target=1, k.min=2)

tclust1$cluster

#--------------------------------------- Tsne
library(tsne)
library(plotly)

set.seed(0)
tsne <- tsne(t(vv), initial_dims = 3, max_iter = 100)
tsne <- data.frame(tsne)
rownames(tsne) <- rownames(vv)
plot(tsne$X1, tsne$X2)
pdb <- cbind(tsne,iris$Species)
options(warn = -1)
fig <-  plot_ly(data = tsne ,x =  ~X1, y = ~X2, type = 'scatter', mode = 'markers')

fig <- fig %>%
  layout(
    plot_bgcolor = "#e5ecf6"
  )

fig
