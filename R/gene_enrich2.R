

enrich_genes <- function(gene_df = res4,
                         gene_col,
                         category = "KEGG",
                         organism = "rat",
                         path_n = 10,
                         score_threshold = 200,
                         version = "11.5") {
  string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold, version = version)
  gene_map <- string_db$map(as.data.frame(gene_df), {{gene_col}}, removeUnmappedRows = TRUE)
  enrichment <- string_db$get_enrichment(gene_map, category = category)
  total_path <- nrow(enrichment)
  if (path_n > total_path) {
    cli_alert_warning(c("You have selected {path_n} pathway{?s}, but only enriched {total_path} pathway{?s}."))
    kegg_res <- enrichment[1:total_path,]
  } else{
  kegg_res <- enrichment[1:path_n,]
  }
  path_tab <- kegg_res %>%
    dplyr::select(c(category, term, description, number_of_genes, p_value, fdr)) %>%
    dplyr::rename(ID= term, path_name = description, n_genes = number_of_genes) %>%
    dplyr::as_tibble()
  if (total_path > 0) {
    gene_list <- strsplit(kegg_res$preferredNames, ",")
    genes_path <- unlist(lapply(gene_list, function(x) {paste(x, collapse = ", ")}))
    path_tab <- path_tab %>%
      dplyr::mutate(genes = genes_path)
  } else {
    path_tab <- path_tab %>%
      dplyr::mutate(genes = kegg_res$preferredNames)
  }
  return(path_tab)
}

get_enrich <- function(gene_df,
                       gene_col,
                       pvalue_col,
                       category = "KEGG",
                       organism = "rat",
                       path_n = 5,
                       score_threshold = 200,
                       version = "11.5") {

  enrc_tab <- enrich_genes(
    gene_df = gene_df,
    gene_col = gene_col,
    category = category,
    organism = organism,
    path_n = path_n,
    score_threshold = score_threshold,
    version = version
  )
  gene_list <- strsplit(enrc_tab$genes, split = "\\, ")
  gene_n <- as.list(enrc_tab$n_genes)
  rep_path <- mapply(function(x, y) {rep(y, x)},x = gene_n, y = enrc_tab$path_name)
  gn_col <- gene_df %>%
    dplyr::select({{gene_col}})
  gene_name <- as.vector(unlist(gn_col))

  edge_data <- tibble(from = unlist(rep_path), to = unlist(gene_list), weight = NA) %>%
    filter(to %in% gene_name)
  enrc_genes <- unique(edge_data$to)
  if(all(grepl(".*e-", enrc_tab$p_value))){
    path_size <- ceiling(scales::rescale(as.numeric(gsub(".*e-","", enrc_tab$p_value)), to = c(10, 15)))
  } else {
    path_size <- ceiling(scales::rescale(enrc_tab$p_value, to = c(10, 15)))
  }
  path_data <- enrc_tab %>%
    dplyr::select(c(path_name, n_genes, p_value)) %>%
    dplyr::mutate(nodes = "pathway", size = path_size) %>%
    dplyr::rename(name = path_name, count = n_genes)
  gene_df2 <- gene_df %>%
    dplyr::filter(gene_name %in% enrc_genes)
  p_col <- gene_df2 %>%
    dplyr::select({{pvalue_col}})
  pvalues <- as.vector(unlist(p_col))
  if(all(grepl(".*e-", pvalues))){
    gene_size <- ceiling(scales::rescale(as.numeric(gsub(".*e-","", pvalues)), to = c(3, 7)))
  } else {
    gene_size <- ceiling(scales::rescale(pvalues, to = c(3, 7)))
  }
  gene_pval <- gene_df2 %>%
    dplyr::select(c({{gene_col}},{{pvalue_col}})) %>%
    dplyr::mutate(nodes = "gene", size = gene_size) %>%
    dplyr::rename(name = {{gene_col}}, p_value = {{pvalue_col}})
  gene_count <- edge_data %>%
    dplyr::rename(name = to) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(count = n())

  gene_data <- left_join(gene_count, gene_pval, by = "name")
  vertx_data <- bind_rows(path_data, gene_data)


  net_data <- list(vertices = vertx_data, edges = edge_data)
  N <- sum(net_data$vertices$nodes == "gene")
  for(i in 1:nrow(net_data$edges)) {
    n1 <- net_data$vertices$count[net_data$vertices$name %in% net_data$edges[i, ]$from]
    n2 <- net_data$vertices$count[net_data$vertices$name %in% net_data$edges[i, ]$to]
    net_data$edges$weight[i] <- N / (n1 * n2)
  }
  return(net_data)
}


cluster_enrich <- function(gene_df,
                           net_data,
                           cluster = 2,
                           gene_col = "gene_name",
                           pvalue_col = "group",
                           category = "KEGG",
                           organism = "rat",
                           path_n = 10,
                           score_threshold = 200,
                           version = "11.5") {
  cluster_df <- net_data$vertices %>%
    dplyr::filter(community == cluster)
  gcluster_df <- gene_df %>%
    dplyr::filter(gene_name %in% cluster_df$name)
  enrc_data <- get_enrich(gene_df = gcluster_df,
                          gene_col = gene_col,
                          pvalue_col = pvalue_col,
                          category = category,
                          organism = organism,
                          path_n =path_n,
                          score_threshold = score_threshold,
                          version = version)
  return(enrc_data)
}

net_enrich <- function(gene_df = res4,
                       net_df = net_edgebet$vertices,
                       gene_col = "gene_name",
                       category = "KEGG",
                       organism = "rat",
                       path_n = 10,
                       score_threshold = 200,
                       version = "11.5") {

  cluster_data <- net_df %>%
    select(name, gene_class) %>%
    group_by(gene_class) %>%
    do(gene_tab = gene_df[gene_df$gene_name %in% .$name,])
  gene_net <- cluster_data$gene_tab
  names(gene_net) <- cluster_data$gene_class
  enrich_all <- lapply(gene_net, function(x) {
    enrich_genes(
      gene_df = x,
      gene_col = gene_col,
      category = category,
      organism = organism,
      path_n = path_n,
      score_threshold = score_threshold,
      version = version
    )
    })
  return(enrich_all)
}

compare_clenrich <- function(gene_df = res4,
                             net_df = net_edgebet$vertices,
                             gene_col = "gene_name",
                             category = "KEGG",
                             organism = "rat",
                             path_n = 10,
                             score_threshold = 200,
                             version = "11.5") {
  full_res <- enrich_genes(gene_df = gene_df,
                           gene_col = gene_col,
                           category = category,
                           organism = organism,
                           path_n = path_n,
                           score_threshold = score_threshold,
                           version = version)
  cluster_res <- net_enrich(gene_df = gene_df,
                            net_df = net_df,
                            gene_col = gene_col,
                            category = category,
                            organism = organism,
                            path_n = path_n,
                            score_threshold = score_threshold,
                            version = version)
  cluster_path <- lapply(cluster_res, function(x) x$path_name)
  enrich_pct <- lapply(cluster_path, function(x) sum(full_res$path_name %in% x)/path_n)
  enrich_res <- data.frame(Cluster = names(enrich_pct), value = unlist(enrich_pct))
  return(enrich_res)
}

compare_netenrich <- function(...,
                              gene_df = res4,
                              gene_col = "gene_name",
                              category = "KEGG",
                              organism = "rat",
                              path_n = 10,
                              score_threshold = 200,
                              version = "11.5") {
  net_list <- list2(...)
  if (length(net_list) == 1 && is_bare_list(net_list[[1]])) {
    net_list <- net_list[[1]]
  }
  if (is.null(names(net_list))) {
    names(net_list) <- paste("Method", LETTERS[1:length(net_list)], sep = "_")
  }
  enrich_all <- lapply(net_list, function(x) {
    compare_clenrich(gene_df = gene_df,
                     net_df = x,
                     gene_col = gene_col,
                     category = category,
                     organism = organism,
                     path_n = path_n,
                     score_threshold = score_threshold,
                     version = version)})
  enrich_all <- mapply(function(x, y) x %>% mutate(Method = y),
         enrich_all, names(enrich_all), SIMPLIFY = FALSE)
  enrich_res <- do.call(rbind, enrich_all)
  return(enrich_res)
}
#---------------------  function for enrich plot ----------------------
plot_enrich <- function(enrc_data, plot_layout = "fr",
                        title = NULL, title_size = rel(1),
                        tag = NULL, tag_size =rel(1)) {
  enrc_net <- igraph::graph_from_data_frame(
    d = enrc_data$edges,
    vertices = enrc_data$vertices,
    directed = TRUE
  )
  col_ver <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF9F00", "#FF00FF",
               "#8B4513", "#556B2F", "#6A3D9A", "yellow3", "paleturquoise4")
  col_edg <- c("deepskyblue", "lightgreen", "#FB9A99", "#FDBF6F", "#DDA0DD",
               "#CD853F", "#CAFF70", "#CAB2D6", "#FFFF99", "paleturquoise3")
  n_path <- sum(enrc_data$vertices$nodes == "pathway")
  n_gene <- sum(enrc_data$vertices$nodes == "gene")

  ggcol <- colorRampPalette(c("azure1", "azure4"))
  logp <- ceiling(-log10(igraph::V(enrc_net)$p_value))[igraph::V(enrc_net)$nodes == "gene"]
  maxp <- plyr::round_any(max(logp), 5, f = ceiling)
  minp <- plyr::round_any(min(logp), -5, f = ceiling)
  gene_col <- ggcol(maxp-minp)[logp-minp]
  igraph::V(enrc_net)$color <- c(col_ver[1:n_path], gene_col)
  igraph::E(enrc_net)$width <- floor(scales::rescale(igraph::E(enrc_net)$weight, to = c(1, 3)))
  igraph::E(enrc_net)$color <- col_edg[match(enrc_data$edges$from, igraph::V(enrc_net)$name)]
  pp <- ggraph::ggraph(enrc_net, layout = plot_layout) +
    ggraph::geom_edge_arc(aes(colour = color, width= 1.5),strength = 0.2, alpha = 0.5) +
    ggraph::geom_node_point(size = c(igraph::V(enrc_net)$size[1:n_path], rep(4, n_gene)),
                            color = igraph::V(enrc_net)$color) +
    ggraph::geom_node_text(aes(label = c(rep("", n_path),attr(igraph::V(enrc_net), "name")[(n_path +1):(n_path + n_gene)])), size = 3, repel = TRUE) +
    labs(title = title,tag = tag) +
    theme(plot.title = element_text(face = "bold", size = rel(title_size),hjust = 0.5),
          plot.tag = element_text(face = "bold", size = tag_size),
          panel.background = element_blank(),
          plot.background = element_rect(colour = NA),
          legend.position = "none"
    )

  path_name <- factor(igraph::V(enrc_net)$name[enrc_data$vertices$nodes == "pathway"],
                      levels= igraph::V(enrc_net)$name[enrc_data$vertices$nodes == "pathway"])
  path_size <- igraph::V(enrc_net)$size[enrc_data$vertices$nodes == "pathway"]
  path_col <- igraph::V(enrc_net)$color[enrc_data$vertices$nodes == "pathway"]
  df_leg <- data.frame(x=1:n_path,y=1:n_path, path_name= path_name)
  leg_plot <- ggplot(df_leg,aes(x=x,y=y,color=path_name, size = path_name))+
    geom_point()+
    scale_color_manual(values=path_col)+
    scale_size_manual(values =path_size/2)+
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size= 12),
      legend.margin = margin(0,unit = "mm"),
      legend.box.margin= margin(0,unit = "mm"),
      legend.spacing.y = unit(-0.3, 'mm'),
      legend.spacing.x = unit(3, 'mm'),
      legend.key=element_blank())+
    guides(color = guide_legend(ncol = ifelse(n_path > 5, 2, 1), byrow = TRUE))
  leg_path <- common_legend(leg_plot)

  df_gene <- data.frame(x=logp,y=1:length(logp))
  leg_gplot <- ggplot(df_gene,aes(x=x,y=y,color=x))+
    geom_point()+
    labs(color = "")+
    scale_colour_gradient(low = "azure1", high = "azure4")
  leg_gene <-common_legend(leg_gplot)

  path_logp <- ceiling(-log(igraph::V(enrc_net)$p_value[enrc_data$vertices$nodes == "pathway"]))
  df_psize <- data.frame(
    x= ceiling(c(min(path_logp), median(path_logp), max(path_logp))),
    y= 1:3,
    z= ceiling(c(min(path_size/2), median(path_size/2), max(path_size/2))))
  leg_ppath <- ggplot(df_psize,aes(x=x,y=y,size=z), color = "gray50")+
    geom_point(aes(size=x))+
    labs(size = "")+
    theme(
      legend.text = element_text(size= 10),
      legend.margin = margin(0,unit = "mm"),
      legend.box.margin= margin(0,unit = "mm"),
      legend.spacing.y = unit(-0.3, 'mm'),
      legend.spacing.x = unit(1, 'mm'),
      legend.key=element_blank())
  leg_psize <-common_legend(leg_ppath)

  ppp <- grid.arrange(pp,
                      gridExtra::arrangeGrob(leg_path,
                      gridExtra::arrangeGrob(grobs = list(leg_psize, leg_gene),
                                             widths = c(5,5), ncol = 2,
                                             top = textGrob(expression("-log"[10]*" (p-value)"))),
                      ncol=2,widths=c(10, 2)),
                      nrow=2,heights=c(10, 3.5))

  return(ppp)
}
