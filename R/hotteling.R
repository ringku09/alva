data_matrix <-naov_data2$compound
# Perform PCA
pca_result <- prcomp(data_matrix, center = TRUE,scale. = TRUE)
# Calculate radius
confidence_level <- 0.95
F <- qf(confidence_level, 2, nrow(data_matrix) - 1)
radius <- sqrt(F * (nrow(data_matrix) - 1) * ncol(data_matrix) / (nrow(data_matrix) - ncol(data_matrix)))

# Generate circle
theta <- seq(0, 2 * pi, length.out = 200)
circle_x <- radius * cos(theta)
circle_y <- radius * sin(theta)
circle_data <- data.frame(x = circle_x, y = circle_y)

# Create a new data frame for PCA results
df_pca <- data.frame(pca_result$x)

# Calculate distance from origin for each point
df_pca$distance <- sqrt(df_pca$PC1^2 + df_pca$PC2^2)

# Create a new variable to indicate if point is inside or outside the circle
df_pca$inside_circle <- ifelse(df_pca$distance <= radius, "EE", "DE")

# Generate plot
p <- ggplot(df_pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = inside_circle), size = 2) +
  geom_path(data = circle_data, aes(x = x, y = y), colour = "red", linetype = "solid", linewidth =1) +
  labs(title = "",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(name = "", values = c("EE" = "white", "DD" = "gray10")) +
  theme_Publication()+theme(legend.position="none")


# Print the plot
ggsave(filename = "E:/lavana/results/Figure3.jpeg", plot = p, width = 190, height = 150, dpi = 300, units = "mm")

table(df_pca$inside_circle)
