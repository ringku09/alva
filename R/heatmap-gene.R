

#' Heatmap of average gene expresion data
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' `plot_heatmap()` will plot heatmap of training data.
#'
#'
#' @param group_var
#' @inheritParams plot_gene
#'
#' @return
#' A plot
#' @export
#'
#' @examples
plot_heatmap <- function(...,
                         expr_data,
                         attr_data,
                         probe_id,
                         group_var = "dose",
                         organism = "rat",
                         multicore = FALSE,
                         store = FALSE,
                         output_dir = missing_arg(),
                         error_call = caller_env()) {
    gene <- probes2genes(probe_id, organism = organism)
    gene_name <- block_fst(gene$GENENAME)
    avg_data <- get_avgFC(...,
                          expr_data = expr_data,
                          attr_data = attr_data,
                          probe_id = probe_id,
                          multicore = multicore,
                          store = store,
                          output_dir = output_dir,
                          error_call = error_call)
    # Cluster compounds based on average expression
    avg_mat <- tidyr::pivot_wider(avg_data$time_mlfc %>% dplyr::select(-group),
                                  names_from = compound, values_from = avg_x)
    avg_mat <- as.matrix(avg_mat%>% dplyr::select(-c(dose,time)))
    clust <- hclust(dist(t(avg_mat)))
    avg_data$time_mlfc$compound <- factor(avg_data$time_mlfc$compound,
                                          levels=colnames(m)[clust$order])
    if (identical(group_var, "dose")) {
      p <- ggplot(avg_data$time_mlfc, aes(x = time, y = compound, fill = avg_x)) +
        geom_tile() +
        facet_grid(group~dose, scales = "free_y", space = "free_y") +
        scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red",
                             space ="Lab", name = expression(bar("log"[2]*" (FC)")))+
        labs(title = glue::glue("{gene_name} ({gene$SYMBOL})"))+
        theme(plot.title = element_text(color="#556B2F", size=15, face="bold", hjust = 0.5))
    } else if (identical(group_var, "time")) {
      p <- ggplot(avg_data$time_mlfc, aes(x = dose, y = compound, fill = avg_x)) +
        geom_tile() +
        facet_grid(group~time, scales = "free_y", space = "free_y") +
        scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red",
                             space ="Lab", name = expression(bar("log"[2]*" (FC):")))+
        labs(title = glue::glue("{gene_name} ({gene$SYMBOL})"))+
      theme(plot.title = element_text(color="#556B2F", size=15, face="bold", hjust = 0.5))
    } else {
      cli_abort(c("Wrong {.arg group_var}.",
                  "x" = "The variable {style_bold(col_red(backtick(group_var)))} is not supported.",
                  "i" = "Please select either `dose` or `time`.")
                , call = error_call)
      }
    p
}


gene_heatmap <- function(...,
                         expr_data,
                         attr_data,
                         probe_id,
                         organism = "rat",
                         com_annotate = TRUE,
                         multicore = FALSE,
                         store = FALSE,
                         tag = NULL,
                         tag_size = rel(1),
                         output_dir = missing_arg(),
                         error_call = caller_env()) {
  gene <- probes2genes(probe_id, organism = organism)
  gene_name <- block_fst(gene$GENENAME)
  avg_data <- get_avgFC(...,
                        expr_data = expr_data,
                        attr_data = attr_data,
                        probe_id = probe_id,
                        multicore = multicore,
                        store = store,
                        output_dir = output_dir,
                        error_call = error_call)
  if (com_annotate) {
    # Cluster compounds based on average expression
    avg_mat <- tidyr::pivot_wider(avg_data$time_mlfc %>% dplyr::select(-group),
                                  names_from = compound, values_from = avg_x)
    avg_mat <- as.matrix(avg_mat%>% dplyr::select(-c(dose,time)))
    clust <- hclust(dist(t(avg_mat)))
    avg_data$time_mlfc$compound <- factor(avg_data$time_mlfc$compound,
                                          levels=colnames(avg_mat)[clust$order])
  }
    p1 <- ggplot(avg_data$time_mlfc, aes(x = time, y = compound, fill = avg_x)) +
      geom_tile() +
      facet_grid(group~dose, scales = "free_y", space = "free_y") +
      scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red",
                           space ="Lab", name = expression(bar("log"[2]*" (FC)")))+
      labs(title = "", tag = tag)+
      theme_Publication(title_size = rel(0.7), axix_size = rel(0.5), lab_size = rel(0.5),
                        x_lab = element_blank(),y_lab = element_blank(),tag_size = tag_size,
                        top_mar = 2)+
      theme(legend.position = "right",
            legend.direction = "vertical")
    p2 <- ggplot(avg_data$time_mlfc, aes(x = dose, y = compound, fill = avg_x)) +
      geom_tile() +
      facet_grid(group~time, scales = "free_y", space = "free_y") +
      scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red",
                           space ="Lab", name = expression(bar("log"[2]*" (FC)")))+
      theme_Publication(title_size = rel(0.7), axix_size = rel(0.5), lab_size = rel(0.5),
                        x_lab = element_blank(),y_lab = element_blank(),
                        top_mar = 2)+
      labs(title = "")+
      theme(legend.position = "right",
            legend.direction = "vertical")
    gene_leg <- com_legend(p2)
    ggg <- gridExtra::grid.arrange(p1 + theme(legend.position = "none"),
                                   p2 + theme(legend.position = "none"),
                                   nrow = 2)
    gene_heat <- gridExtra::grid.arrange(ggg,
                                           gene_leg,ncol=2, widths = c(9,1),
                                           top = grid::textGrob(glue::glue("{gene_name} ({gene$SYMBOL})"),
                                                                gp=grid::gpar(fontsize=16, font = 11)))

    gene_heat
}


tgx_heatmap <- function(...,
                        expr_data,
                        attr_data,
                        probes,
                        dose = "High",
                        time = "24 hr",
                        annotation_compound = FALSE,
                        annotation_legend = TRUE) {
  probes <- unique(probes)
  comps_gr <- list2(...)
  if (length(comps_gr) == 1 && is_bare_list(comps_gr[[1]])) {
    comps_gr <- comps_gr[[1]]
  }
  if (is.null(names(comps_gr))) {
    names(comps_gr) <- paste("Group", LETTERS[1:length(comps_gr)], sep = "_")
  }
  space_data <- avg_space(
    comps_gr,
    expr_data = expr_data,
    attr_data = attr_data,
    probes = probes,
    dose = dose,
    time = time)
  sample_expr <- space_data$expr_data
  # set genename
  genes <- probes2genes(probes) %>%
    mutate(SYMBOL = coalesce(SYMBOL, PROBEID))
  rownames(sample_expr) <-  genes$SYMBOL
  sample_expr <- sample_expr[!duplicated(rownames(sample_expr)), ]

  anno_sample <- space_data$attr_data
  tgx_class <- vector(length = nrow(anno_sample))
  for(i in 1:length(comps_gr)) {
    tgx_class[anno_sample$compound_abbr %in% comps_gr[[i]]] <- names(comps_gr)[i]
  }
  tgx_class <- factor(tgx_class, levels = names(comps_gr))
  anno_sample <- anno_sample %>%
    dplyr::mutate(Group = tgx_class, .before = compound_abbr) %>%
    dplyr::rename(Compound = compound_abbr) %>%
    select(Group, Compound)
  if(annotation_compound) {
    anno_attr <- anno_sample %>%
      select(Group, Compound)
  } else {
    anno_attr <- anno_sample %>%
      select(Group)
  }
  anno_attr <- as.data.frame(anno_attr)
  rownames(anno_attr) <- colnames(sample_expr)
  col_group <- gr_cols <- c("firebrick3", "springgreen", "dodgerblue3", "orange3")



  anno_gene <- up_down(comps_gr,
                       expr_data= expr_data,
                       attr_data = attr_data,
                       probes = probes)
  anno_gene <- anno_gene %>%
    dplyr::mutate(gene_id = genes$SYMBOL) %>%
    dplyr::rename(Genes = Expression) %>%
    dplyr::select(gene_id, Genes) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
  anno_gene <- as.data.frame(anno_gene %>% dplyr::select(Genes))
  rownames(anno_gene) <- rownames(sample_expr)
  anno_colour = list(
    Group = col_group[1 : length(levels(anno_attr$Group))],
    Genes = c('yellowgreen', 'pink')
  )
  names(anno_colour$Group) = levels(anno_attr$Group)
  names(anno_colour$Genes) = levels(anno_gene$Genes)

  # genes <- probes2genes(gaogene) %>%
  #   mutate(SYMBOL = coalesce(SYMBOL, PROBEID))
  # zz <- zz[!duplicated(rownames(zz)), ]
  # rownames(zz) <-  genes$SYMBOL
  #

  p <- pheatmap(sample_expr, scale = "row", color = hcl.colors(1000, "PiYG"),
                annotation_colors = anno_colour,
                annotation_row = anno_gene,
                annotation_col = anno_attr,
                annotation_legend = annotation_legend)
  return(p)
}
