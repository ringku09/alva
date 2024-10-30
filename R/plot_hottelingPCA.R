get_pcaT2 <- function(pca_data, gene_tab, k=2, alpha = 0.05, pts = 200) {
  pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  contb <- round(pca$sdev^2/sum(pca$sdev^2)*100, 1)
  pca_scores <- as.matrix(pca$x)
  sigma2 <- pca$sdev^2
  n <- nrow(pca_scores)
  p <- 1 - alpha
  # 99% and 95% confidence limit
  F_dist <- stats::qf(p, df1 = k, df2 = (n-k))
  tsq_critical <- (k*(n-1)/(n-k)) * F_dist
  # Hotelling ellipse semi-axes
  a <- sqrt(tsq_critical*sigma2[1])
  b <- sqrt(tsq_critical*sigma2[2])
  tsq <- (n-k)/(k*(n-k)) * diag(pca_scores[,1:k] %*% diag(1/sigma2[1:k]) %*% t(pca_scores[,1:k]))
  F_stat <- ((n - 1) * k / (n - k)) * tsq / n
  p_value <- 1 - pf(tsq, k, n-k)
  log10p = -log10(p_value)
  # de_cutoff <- (pca_scores[,1]^2 / a^2 + pca_scores[,2]^2 / b^2) <= 1
  ee_de <- ifelse(tsq < F_dist, "EE", "DE")
  gene_symbol <- gene_tab$gene_name[match(rownames(pca_scores), gene_tab$probe_id)]
  tox_log10p <- gene_tab$log10p[match(rownames(pca_scores), gene_tab$probe_id)]
  gene_name <- gene_tab$gene_fnane [match(rownames(pca_scores), gene_tab$probe_id)]
  pca_data <- tibble::tibble(probe_id = rownames(pca_scores), gene_symbol = gene_symbol,
                             gene_name = gene_name, PC1 = pca_scores[,1], PC2 = pca_scores[,2],
                             T2 = tsq, p_value = p_value,log10p =log10p, gene = ee_de)
  # Coordinate points
  pp <- seq(0, 2 * pi, length.out = pts)
  ellips_x <- a*cos(pp)
  ellips_y <- b*sin(pp)
  ellips_data <- tibble::tibble(x = ellips_x, y = ellips_y)
  return(list(pca_data = pca_data, ellips_data = ellips_data, T2_critical = tsq_critical, pc_contb = contb))
}



plot_pcaellips <- function(pca_data, gene_tab, alpha = 0.05,top = 30, pts = 200) {
  pca_ellpsres <- get_pcaT2(pca_data = pca_data,gene_tab = gene_tab, k=2, alpha = alpha, pts = pts)
  pca_elps <- pca_ellpsres$pca_data
  top_idx <- order(pca_elps$T2, decreasing = TRUE)[1:top]
  pca_elps$lab <- ""
  pca_elps$lab[top_idx] <- pca_elps$gene_symbol[top_idx]
  p <- ggplot(pca_elps, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = gene), size = 2) +
    ggrepel::geom_text_repel(aes(x = PC1, y = PC2,label=lab)) +
    geom_path(data = pca_ellpsres$ellips_data, aes(x = x, y = y), colour = "red",
              linetype = "solid", linewidth =1,inherit.aes = FALSE) +
    labs(title = "",
         x = glue("PC1 ({pca_ellpsres$pc_contb[1]}%)"),
         y = glue("PC2 ({pca_ellpsres$pc_contb[2]}%)")) +
    scale_color_manual(name = "", values = c("EE" = "gray90", "DD" = "gray10")) +
    theme_Publication()+theme(legend.position="none")
  return(p)
}



plot_pcamhat <- function(pca_data, gene_tab, k = 2, alpha = 0.05) {
  pca_ellpsres <- get_pcaT2(pca_data = pca_data, gene_tab = gene_tab, k=k, alpha = alpha)
  pca_elps <- pca_ellpsres$pca_data
  p <- pca_elps %>%
    mutate(obs = 1:nrow(pca_elps)) %>%
    ggplot(aes(x = obs, y = T2)) +
    #geom_point(aes(fill = T2), shape = 21, size = 3, color = "gray90") +
    geom_path(aes(x = obs, y = T2), size = .5, colour = "gray50") +
    #scale_fill_gradient(low = "white", high = "red", guide = FALSE) +
    geom_hline(yintercept = pca_ellpsres$T2_critical, linetype = "dashed", color = "darkred", size = .5) +
    annotate(geom = 'text', label = paste('  # Selected genes', "=", sum(pca_elps$gene == "DE")),
             x = -Inf, y = Inf, hjust = 0, vjust = 1)+
    # geom_text(x = nrow(tt$pca_data), y = max(tt$pca_data$T2), label = paste('genes', "==", ngene),
    #           parse = TRUE, vjust = -1, color = "darkred")+
    labs(x = expression("Gene (ordered by -log"[10]*"P-value)"), y = expression("Hotelling’s T"^2),
         fill = expression("-log"[10]*" (P-value)")) +
    theme_Publication()
  return(p)
}


