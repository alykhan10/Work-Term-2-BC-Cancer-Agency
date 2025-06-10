# Install required packages
install.packages("pamr")
install.packages("tidyverse")

# Load libraries
library(pamr)
library(tidyverse)

# Load gene expression data
expr_df <- read_csv("BLCA.csv")

# Convert to matrix with gene names as rownames
expr_mat <- expr_df %>%
  column_to_rownames(var = colnames(expr_df)[1]) %>%
  as.matrix()

# Load sample subtype labels
labels_df <- read_csv("BLCA_labels.csv")

# Reorder expression columns to match label sample IDs
expr_mat <- expr_mat[, labels_df$sample_id]

# Format data for pamr
pamr_data <- list(
  x = expr_mat,
  y = as.factor(labels_df$subtype),
  geneid = rownames(expr_mat),
  genenames = rownames(expr_mat)
)

# Train nearest-centroid classifier
pamr_model <- pamr.train(pamr_data)

# Cross-validate to find optimal threshold
cv_results <- pamr.cv(pamr_model, pamr_data)
pamr.plotcv(cv_results)

# Predict subtypes using a chosen threshold
threshold <- 0.5
predicted_subtypes <- pamr.predict(pamr_model, pamr_data$x, threshold = threshold)

# Confusion matrix to evaluate performance
conf_matrix <- table(Predicted = predicted_subtypes, Actual = pamr_data$y)
print(conf_matrix)

# List top genes contributing to subtype separation
top_genes <- pamr.listgenes(pamr_model, pamr_data, threshold = threshold)
head(top_genes)