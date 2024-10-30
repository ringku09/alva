# Required Libraries
library(ggplot2)
library(plotly)
library(ggfortify)
library(ellipse)

# Generate synthetic data
data <- matrix(rnorm(100*4, mean=0, sd=1), 100, 4)
colnames(data) <- c("Var1", "Var2", "Var3", "Var4")

# Run PCA
data_matrix <-bb$alpha

pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)

# Extract the scores and get their names
scores <- as.data.frame(pca_result$x)
pc_names <- names(scores)

# Extract loadings
loadings <- pca_result$rotation[, 1:2]

# Calculate covariance matrix and get ellipse
cov_mat <- cov(scores)/500
ell_data <- as.data.frame(ellipse(cov_mat, centre=c(0, 0)))
names(ell_data) <- c("x", "y") # explicitly name the columns

# Plotting using ggplot
p <- ggplot(data = scores, aes(x = pc_names[1], y = pc_names[2])) +
  geom_point(aes(color = "red")) +
  geom_path(data = ell_data, aes(x = x, y = y), color = 'blue') +
  ggtitle("PCA Biplot with Ellipse")

# Add loadings
p + geom_segment(data = as.data.frame(loadings), aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
                   arrow = arrow(type = "closed", length = unit(0.2, "inches")),
                   color = "black", size = 0.5)



library(plotly)

df <- iris[1:4]
pca_res <- prcomp(df, scale. = TRUE)

p <- autoplot(pca_result, data = dt, colour = 'red',
              loadings = TRUE, loadings.colour = 'blue',
              loadings.label = TRUE, loadings.label.size = 3)+
  geom_path(data = ell_data, aes(x = x, y = y), color = 'blue')

ggplotly(p)
####################################################


# Generate example data
data_matrix <-naov_data2$compound
# Perform PCA
pca_result <- prcomp(data_matrix, center = TRUE,scale. = TRUE)

confidence_level <- 0.99
F <- qf(confidence_level, 2, nrow(data_matrix) - 1)
radius <- sqrt(2 * F * (nrow(data_matrix) - 1) * (nrow(data_matrix) - 1) / (nrow(data_matrix) * (nrow(data_matrix) - 2)))

theta <- seq(0, 2 * pi, length.out = 100)
circle_x <- radius * cos(theta)
circle_y <- radius * sin(theta)
circle_data <- data.frame(x = circle_x, y = circle_y)

p <- ggplot(data.frame(pca_result$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = "Score"), size = 2) +
  geom_path(data = circle_data, aes(x = x, y = y), colour = "red", linetype = "dashed") +
  labs(title = "PCA Biplot with Hotelling's T^2 Circle",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(name = "", values = c("Score" = "gray", "Loading" = "red")) +
  theme_minimal()

print(p)

##--------------------------------------------- hotteling

library(ggforce)
pca_mod <- data_matrix %>%
  PCA(scale.unit = FALSE, graph = FALSE, ncp = 9)

pca_scores <- pca_mod %>%
  pluck("ind", "coord") %>%
  as_tibble() %>%
  print()


res_2PCs <- ellipseParam(data = pca_scores, k = 2, pcx = 1, pcy = 2)

a1 <- pluck(res_2PCs, "Ellipse", "a.99pct")
b1 <- pluck(res_2PCs, "Ellipse", "b.99pct")

a2 <- pluck(res_2PCs, "Ellipse", "a.95pct")
b2 <- pluck(res_2PCs, "Ellipse", "b.95pct")

T2 <- pluck(res_2PCs, "Tsquare", "value")

coord_2PCs_99 <- ellipseCoord(data = pca_scores, pcx = 1, pcy = 3, conf.limit = 0.99, pts = 200)
coord_2PCs_95 <- ellipseCoord(data = pca_scores, pcx = 1, pcy = 3, conf.limit = 0.95, pts = 200)
coord_2PCs_90 <- ellipseCoord(data = pca_scores, pcx = 1, pcy = 3, conf.limit = 0.90, pts = 200)


pca_scores %>%
  ggplot(aes(x = Dim.1, y = Dim.2)) +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = .5, linetype = "dotted", fill = "white") +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a2, b = b2, angle = 0), size = .5, linetype = "dashed", fill = "white") +
  geom_point(aes(fill = T2), shape = 21, size = 1, color = "black") +
  scale_fill_viridis_c(option = "viridis") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = .2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = .2) +
  labs(title = "Scatterplot of PCA scores", subtitle = "PC1 vs. PC2", x = "PC1", y = "PC2", fill = "T2", caption = "Figure 1: Hotelling's T2 ellipse obtained\n using the ellipseParam function") +
  theme_grey()


ggplot() +
  geom_ellipse(data = coord_2PCs_99, aes(x0 = x, y0 = y, a = 1, b = 1, angle = 0), size = .5, color = "black", linetype = "solid") +
  geom_ellipse(data = coord_2PCs_95, aes(x0 = x, y0 = y, a = 1, b = 1, angle = 0), size = .5, color = "darkred", linetype = "solid") +
  geom_ellipse(data = coord_2PCs_90, aes(x0 = x, y0 = y, a = 1, b = 1, angle = 0), size = .5, color = "darkblue", linetype = "solid") +
  geom_point(data = pca_scores, aes(x = Dim.1, y = Dim.3, fill = T2), shape = 21, size = 1, color = "black") +
  scale_fill_viridis_c(option = "viridis") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = .2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = .2) +
  labs(title = "Scatterplot of PCA scores", subtitle = "PC1 vs. PC3", x = "PC1", y = "PC3", fill = "T2", caption = "Figure 2: Hotelling's T2 ellipse obtained\n using the ellipseCoord function") +
  theme_bw() +
  theme(panel.grid = element_blank())
