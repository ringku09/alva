
#' Setup STRING database
#'
#' The function `setup_stringdb()` is used to setup STRING database.
#'
#' @param organism Sample organism used for experiment ("human" or "rat")
#' @param score_threshold A threshold for the combined scores of the interactions  (default 200)
#' @param version Version of STRING database. A string, e.g., "11.5" (the default).
#'
#' @return Search Tool for the Retrieval of Interacting Genes/Proteins
#' @export
#'
#' @examples
#' \dontrun{
#' setup_stringdb(organism = "rat", score_threshold = 200)
#' }
setup_stringdb <- function(organism = "rat", score_threshold = 200, version = "11.5") {
  if (identical(organism, "human")) {
    taxa = 9606
  } else if (identical(organism, "rat")) {
    taxa = 10116
  }
  string_db <- STRINGdb::STRINGdb$new(
    version = version,
    species = taxa,
    network_type = "full",
    score_threshold = score_threshold,
    input_directory = ""
    )
  return(string_db)
}


#' Get network data
#'
#' The function `get_netdata()` is used to extract network data from gene table.
#'
#' @param gene_df A data frame of gene information or a `tgxtool` class of object.
#' @param gene_col The name of the column for gene name in `gene_df`.
#' @param p.value_col The name of the column for p-values in `gene_df`.
#' @param cluster_method Cluster algorithm to find community. You can choose between "fastgreedy",
#' "walktrap", "spinglass" and "edge.betweenness" (the default is `NULL`).
#' @inheritParams setup_stringdb
#'
#' @return A list contain two data frames
#'  1) `vertices`: A data frame of information about nodes.
#'  2) `edges`: A data frame of information about edges.
#' @export
#'
#' @examples
#' \dontrun{
#' net_data <- get_netdata(
#'   res3[1:200,],
#'   gene_col = "gene_name",
#'   p.value_col = "group",
#'   organism = "rat",
#'   score_threshold = 400,
#'   cluster_method = "edge.betweenness"
#' )
#' }
get_netdata <- function(gene_df=res4,
                        probe_col = "probe_id",
                        gene_col = "gene_name",
                        p.value_col = "group",
                        cluster_method = NULL,
                        organism = "rat",
                        score_threshold = 200,
                        version = "11.5") {
# test `gene_col` is string  and `p.value_col` is p-values and both contain
# set test for  cluster_method
  string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold, version = version)
  gene_map <- string_db$map(as.data.frame(gene_df), {{gene_col}}, removeUnmappedRows = TRUE)
  gene_map <- gene_map[!duplicated(gene_map$STRING_id),]
  net <- string_db$get_subnetwork(gene_map$STRING_id)
  net_df <- suppressWarnings(igraph::as_data_frame(net, what = "both"))
  net_df$edges <- net_df$edges %>%
    dplyr::mutate(
      from = dplyr::pull(gene_map[{{gene_col}}])[match(net_df$edges$from, gene_map$STRING_id)],
      to = dplyr::pull(gene_map[{{gene_col}}])[match(net_df$edges$to, gene_map$STRING_id)]
      ) %>%
    distinct()
  rm_duplt <- net_df$vertices %>%
    dplyr::mutate(genes = dplyr::pull(gene_map[{{gene_col}}])[match(net_df$vertices$name, gene_map$STRING_id)]) %>%
    distinct(genes, .keep_all = TRUE)
  net_df$vertices <- net_df$vertices %>%
    dplyr::mutate(name = dplyr::pull(gene_map[{{gene_col}}])[match(net_df$vertices$name, gene_map$STRING_id)]) %>%
    distinct()
  network <- igraph::graph_from_data_frame(d = net_df$edges, vertices = net_df$vertices, directed = TRUE)
  deg <- igraph::degree(network)            # Degree centrality
  clo <- igraph::closeness(network)         # Closeness centrality
  bet <- igraph::betweenness(network)       # Betweenness centrality
  eig <- igraph::evcent(network)$vector     # Eigenvector centrality
  p.values <- gene_map %>%
    dplyr::select({{p.value_col}})
  probe <- gene_map %>%
    dplyr::select({{probe_col}})
  net_df$vertices <- net_df$vertices %>%
    dplyr::mutate(
      probes = probe[match(net_df$vertices$name, dplyr::pull(gene_map[{{gene_col}}])), ],
      p.values = p.values[match(net_df$vertices$name, dplyr::pull(gene_map[{{gene_col}}])), ],
      degree = deg,
      betweenness = bet,
      closenes = clo,
      eigenes = eig
      )
  if (!is.null(cluster_method)) {
    cl_list <- string_db$get_clusters(rm_duplt$name, algorithm = cluster_method)
    min_mem <- ceiling(sum(length(gene_map$STRING_id)) * 0.05)
    act_gr <- lapply(cl_list, function(x) {length(x) >= min_mem})
    # inact_gr <- length(cl_list) - sum(unlist(act_gr))
    others <- list(unlist(cl_list[!unlist(act_gr)]))
    names(others) <- "Miscellaneous" #paste("others", inact_gr, sep = "-")
    act_comunity <- cl_list[unlist(act_gr)]
    names(act_comunity) <- paste0("Cluster-", 1:length(act_comunity))
    comunity <- append(act_comunity, others)
    if (length(comunity) > 10) {
      comunity9 <- comunity[1 : 9]
      # others_gr <- length(comunity) - 9
      others <- list(as.vector(unlist(comunity[10 : length(comunity)])))
      names(others) <- "Miscellaneous" # paste("others", others_gr, sep = "-")
      comunity <- append(comunity9, others)
    }

    comunity <- lapply(comunity, function(x) {dplyr::pull(gene_map[{{gene_col}}])[match(x, gene_map$STRING_id)]})
    cl_name <- 1 : length(comunity)
    community <- unlist(lapply(lapply(net_df$vertices$name, function(x) {grep(x, comunity)}), "[",1))
    gene_class <- names(comunity)[unlist(lapply(lapply(net_df$vertices$name, function(x) {grep(x, comunity)}), "[",1))]
    net_df$vertices <- net_df$vertices %>%
      dplyr::mutate(community = community,
                    gene_class = factor(gene_class, levels = names(comunity))
                    )
  }
  return(net_df)
}

#' Extract hub network data
#'
#' The function `get_hubdata()` is used to extract hub network data from full network data.
#'
#' @param net_data A network data of class `tgxtool`.
#' @param condition A character string specify the condition of hub network.
#' (default is `degree >= 20`). There are four metrics `degree`, `betweenness`, `closenes`,
#'  and `eigenes` can be used in the condition with the five comparison operators `<`, `>`, `<=`,
#'  `>=`, and `==`. At the end of the condition you must provide appropriate numeric value of metric
#'  used in the condition.
#'
#' @return A network data
#' @export
#'
#' @examples
get_hubdata <- function(net_data, condition = "degree >= 10", error_call = caller_env()) {
  condition <- gsub(" ", "", condition)
  attrb <- gsub("[^(A-Za-z)]", "", condition)
  value <- as.numeric(gsub("[^0-9.-]", "", condition))
  cond <- gsub("[a-zA-Z0-9.]", "", condition)

  if (!any(sapply(list("degree", "betweenness", "closenes", "eigenes"), FUN = identical, attrb))) {
    cli_abort(c("Metric used in condition must be a valid metric.",
                "x" = "Input {style_bold(col_red(backtick(attrb)))} is not a valid metric.",
                "i" = "Please use either `degree`, `betweenness`, `closenes`, or `eigenes` instread.")
              , call = error_call)
  }
  if (!any(sapply(list(">", "<", ">=", "<=", "=="), FUN = identical, cond))) {
    cli_abort(c("Comparison operator used in condition must be a valid comparison operator.",
                "x" = "Input {style_bold(col_red(backtick(cond)))} is not a valid metric.",
                "i" = "Please use either `<`, `>`, `<=`, `>=`, or `==` instread.")
              , call = error_call)
  }
  if (!is.numeric(value) | is.na(value)) {
    cli_abort(c("The condition must contain numeric value after symbol.",
                "x" = "The numeric value in condition are missing.",
                "i" = "Please provide appropriate numeric value in the condition.")
              , call = error_call)
  }
  valid_cond <- net_data$vertices %>%
    dplyr::select(!!rlang::parse_expr(attrb)) %>%
    unlist()
  if (all(value > valid_cond)) {
    max_val <- max(valid_cond)
    cli_abort(c("The numeric value must be valid with condition.",
                "x" = "The numeric value {value} in condition exceed the real values.",
                "i" = "Please provide value smaller than {max_val} in the condition.")
              , call = error_call)
  }
  net_data$vertices <- net_data$vertices %>%
    dplyr::filter(!!rlang::parse_expr(condition)) %>%
    droplevels()
  hub_from <- (1 : nrow(net_data$edges))[net_data$edges$from %in% net_data$vertices$name]
  hub_to <- (1 : nrow(net_data$edges))[net_data$edges$to %in% net_data$vertices$name]
  hub_node <- intersect(hub_from, hub_to)
  net_data$edges <- net_data$edges[hub_node, ]
  return(net_data)
}

#' plot default string network
#'
#' @inheritParams get_netdata
#'
#' @return A ppi network plot
#' @export
#'
#' @examples
#' \donttest{
#' plot_stringnet(res3[1:50, ],
#'             gene_col = "gene_name",
#'             organism = "rat",
#'             score_threshold = 200)
#' }
plot_stringnet <- function(gene_df,
                           gene_col = "gene_name",
                           organism = "rat",
                           score_threshold = 200) {
  string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold)
  gene_map <- string_db$map(as.data.frame(gene_df), {{gene_col}}, removeUnmappedRows = TRUE)
  gene_map <- gene_map[!duplicated(gene_map$STRING_id),]
  string_db$plot_network(gene_map$STRING_id)
}


#' Plot protein protein interaction network
#'
#' The function `plot_network()` is used to plot protein protein interaction (PPI) network.
#'
#' @param net_data A network data of class `tgxtool`. The network data must be a list which contain
#' two data frames `vertices` and `edges`. However,  the `vertices` must have column `p.values` and
#' `community`.
#' @param plot_layout
#' @param ... Optional arguments passed to [igraph::plot.igraph()].
#'
#' @return
#' @export
#'
#' @examples
plot_network <- function(net_data, plot_layout = "fr", rm_noedge = TRUE,arc = 0.2,
                         tag =  NULL, tag_size = rel(1), show_legend = TRUE)
  {
  net <- igraph::graph_from_data_frame(d = net_data$edges, vertices = net_data$vertices, directed = TRUE)
  if (rm_noedge) {
    net <- igraph::delete.vertices(net, igraph::degree(net)==0)
  }
  col_ver <- c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
               "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4")
  col_edg <- c("deepskyblue", "lightgreen", "#FB9A99","#CAB2D6", "#DDA0DD",
               "#CD853F", "#CAFF70", "#FFFF99", "#FDBF6F", "paleturquoise3")

  igraph::V(net)$color <- col_ver[igraph::V(net)$community]
  if(all(grepl(".*e-", igraph::V(net)$p.values))){
    igraph::V(net)$size <- ceiling(scales::rescale(
      as.numeric(gsub(".*e-","", igraph::V(net)$p.values)), to = c(5, 10)))
  } else {
    igraph::V(net)$size <- ceiling(scales::rescale(1/igraph::V(net)$p.values, to = c(5, 10)))
  }
  # as.numeric(gsub(".*e-","", igraph::V(net)$p.values))
  igraph::E(net)$width <- scales::rescale(igraph::E(net)$combined_score, to = c(0.1, 1))
  igraph::E(net)$color <- col_edg[igraph::V(net)$community][match(net_data$edges$from, igraph::V(net)$name)]
  pp <- ggraph::ggraph(net, layout = plot_layout) +
    ggraph::geom_edge_arc(aes(colour = color, width= width),strength = arc, alpha = 0.5) +
    ggraph::geom_node_point(size = igraph::V(net)$size,
                            color = igraph::V(net)$color) +
    ggraph::geom_node_text(aes(label = block_fst(igraph::V(net)$name)), size = 3, repel = TRUE) +
    labs(tag = tag)+
    theme(#text = element_text(),
      plot.tag = element_text(face = "bold", size = tag_size),
      panel.background = element_blank(),
      plot.background = element_rect(colour = NA),
      legend.position = "none"
    )
  class_name <- factor(igraph::V(net)$gene_class)
  n_net <- length(levels(class_name))
  net_col <- col_ver[sort(unique(igraph::V(net)$community))]
  df_leg <- data.frame(x=1:n_net,y=1:n_net, net_name= levels(class_name))
  leg_plot <- ggplot(df_leg,aes(x=x,y=y,color=net_name), size = rel(1))+
    geom_point()+
    scale_color_manual(values = net_col)+
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size= 10),
      legend.margin = margin(0,unit = "mm"),
      legend.box.margin= margin(0,unit = "mm"),
      legend.spacing.y = unit(-0.3, 'mm'),
      legend.spacing.x = unit(1, 'mm'),
      legend.key=element_blank())+
    guides(color = guide_legend(override.aes = list(size = 4),
                                nrow = dplyr::case_when(n_net <= 3~ 1,
                                                        n_net > 3 & n_net <= 6 ~ 2,
                                                        n_net > 6 ~ 3), byrow = TRUE))
  leg_net <- common_legend(leg_plot)

  temp_df <- data.frame(x=rnorm(10),y=rnorm(10),z=runif(10,0.1,1),
                        g=sample(letters[1:2],10,replace = TRUE))
  leg_edgwt <- ggplot(temp_df) +
    geom_line(
      aes(x = x, y = y, linewidth = z, group = g),
      lineend = "round"
    ) +
    scale_linewidth("Weight: ", range = c(0.1, 2))+
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size= 10),
      legend.margin = margin(0,unit = "mm"),
      legend.box.margin= margin(0,unit = "mm"),
      legend.spacing.y = unit(-0.3, 'mm'),
      legend.spacing.x = unit(1, 'mm'),
      legend.key=element_blank())
  leg_wtsize <- common_legend(leg_edgwt)

  net_logp <- ceiling(-log10(igraph::V(net)$p.values))
  net_size <- igraph::V(net)$size
  df_nsize <- data.frame(
    x= ceiling(c(min(net_logp), median(net_logp), max(net_logp))),
    y= 1:3,
    z= ceiling(c(min(net_size), median(net_size), max(net_size))))
  leg_pnet <- ggplot(df_nsize,aes(x=x,y=y,size=z), color = "gray50")+
    geom_point(aes(size=x))+
    labs(size = expression("-log"[10]*" (p-value): "))+
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size= 10),
      legend.margin = margin(0,unit = "mm"),
      legend.box.margin= margin(0,unit = "mm"),
      legend.spacing.y = unit(-0.3, 'mm'),
      legend.spacing.x = unit(1, 'mm'),
      legend.key=element_blank())
  leg_nsize <-common_legend(leg_pnet)
   if (show_legend) {
     ppp <- grid.arrange(pp,
                         gridExtra::arrangeGrob(leg_net,
                                                gridExtra::arrangeGrob(grobs = list(leg_nsize, leg_wtsize),
                                                                       heights = c(5,5), nrow = 2),
                                                ncol=2,widths=c(6, 4)),
                         nrow=2,heights=c(10, 1.5))
   } else {
     ppp <- pp
   }
  return(ppp)
}

# legend shape number fix korte hobe

