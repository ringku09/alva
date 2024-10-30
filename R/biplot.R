bi_plot <- function(comps_gr,pca_data) {
  pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  com_expr <- unlist(lapply(str_split(colnames(pca_data), " "),"[",1))
  com_group <- names(comps_gr)[get_class(comps_gr, com_expr)]
  pca_biplot <- fviz_pca_biplot(pca,
                                # Fill individuals by groups
                                geom.ind = "point",
                                arrowsize = 0.8,
                                pointshape = 21,
                                pointsize = 2,
                                fill.ind = "gray90",
                                col.var = com_group,
                                col.ind = "gray90",
                                # Color variable by groups
                                #col.var = rep(c("red","green","blue"),times=c(7,3,7)),
                                legend.title = list(fill = " ", color = "Group: "),
                                repel = FALSE,       # Avoid label overplotting,
                                title = ""
  )+
    scale_colour_Publication()+
    theme_Publication()
  # theme(legend.position="none")
  return(pca_biplot)
}


bi_plot2 <- function(comps_gr,pca_data, dose = "H") {
  pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  com <- unlist(lapply(str_split(colnames(pca_data), " "),"[",1))
  dose_expr <- unlist(lapply(str_split(colnames(pca_data), " "),"[",2))
  doses <- gsub(".*?\\(([A-Za-z]).*", "\\1", dose_expr)
  com_group <- names(comps_gr)[get_class(comps_gr, com)]
  dose_name <- if(identical(dose, "H")) {"High"
    }else if(identical(dose, "M")) {"Middle"
      }else if(identical(dose, "L")) "Low"
  pca_biplot <- fviz_pca_var(pca,
                             arrowsize = 0.8,
                             col.var = com_group,
                             select.var = list(name = colnames(pca_data)[doses %in% dose]),
                             # Color variable by groups
                             #col.var = rep(c("red","green","blue"),times=c(7,3,7)),
                             #legend.title = list(fill = " ", color = "Group: "),
                             repel = FALSE,       # Avoid label overplotting,
                             title = glue("Dose = {dose_name}")
  )+
    scale_colour_Publication()+
    theme_Publication()
  # theme(legend.position="none")
  return(pca_biplot)
}
