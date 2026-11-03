# Load necessary libraries
library(stats)  # For PCA
library(ggplot2) # for plotting

# Assuming 'data' is already loaded. Rows are trials/time, columns are neurons.

# --- PCA on Covariance Matrix -------------------------
covmat <- cov(data)  # Obtain covariance matrix of data
eig_vals <- eigen(covmat)$values  # Eigenvalues
eig_vecs <- eigen(covmat)$vectors  # Eigenvectors

explained_var <- eig_vals / sum(eig_vals) * 100  # Variance explained as percentage

# Perform PCA directly on the data matrix
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)  # Perform PCA
score_PCs <- pca_result$x  # PC scores (principal components)
coef_load <- pca_result$rotation  # Coefficients (loadings)
variancePC <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100  # Variance explained in percentage

# Combine neuron index with the loadings
coef_load <- data.frame(Neuron = 1:ncol(data), coef_load)

# --- Select component to plot -------------------------
comp <- 1  # Choose component to plot

# --- Sort neurons based on their loadings for the selected component
coef_load <- coef_load[order(-coef_load[, comp + 1]), ]

# --- Identify neurons with highest positive and negative loadings
pos_loaded_nrns <- coef_load$Neuron[1:5]  # Top 5 neurons (positive loadings)
mean_posloaded_nrns <- rowMeans(data[, pos_loaded_nrns])  # Mean activity of top 5 positive neurons

neg_loaded_nrns <- coef_load$Neuron[(nrow(coef_load)-4):nrow(coef_load)]  # Top 5 neurons (negative loadings)
mean_negloaded_nrns <- rowMeans(data[, neg_loaded_nrns])  # Mean activity of top 5 negative neurons

# --- Normalize and plot the results -------------------
fplot <- function(score_comp, mean_posloaded_nrns, mean_negloaded_nrns, comp, dLP_idx, iLP_idx) {
  # Normalize components and neuron means
  plot_comp <- scale(score_comp[, comp], center = TRUE, scale = TRUE)
  pos_load <- scale(mean_posloaded_nrns, center = TRUE, scale = TRUE)
  neg_load <- scale(mean_negloaded_nrns, center = TRUE, scale = TRUE)
  
  # Plot PC and average neuron activity
  par(mfrow=c(2, 1))  # Set plotting layout to 2 rows, 1 column
  
  # Plot PC1 and mean spike count of top positively loaded neurons
  plot(plot_comp, type="l", col="blue", ylab="PC1", main="PC1 and Mean Spike Count of Top Positively Loaded Neurons")
  lines(pos_load, col="green")
  
  # Plot PC1 and mean spike count of top negatively loaded neurons
  plot(plot_comp, type="l", col="blue", ylab="PC1", main="PC1 and Mean Spike Count of Top Negatively Loaded Neurons")
  lines(neg_load, col="red")
  
  # 3D plot for PCA components (optional in R, typically would use scatterplot3d or rgl library)
  # You can use libraries like 'scatterplot3d' or 'rgl' for 3D visualization.
  
  # Example with scatterplot3d
  library(scatterplot3d)
  scatterplot3d(score_comp[dLP_idx, 1], score_comp[dLP_idx, 2], score_comp[dLP_idx, 3], color="blue", pch=16)
  scatterplot3d(score_comp[iLP_idx, 1], score_comp[iLP_idx, 2], score_comp[iLP_idx, 3], color="red", pch=16, add=TRUE)
  title("Top 3 Components across Trials\nDelayed LPs (blue) and Immediate LPs (red)")
}

# Call the fplot function
fplot(score_PCs, mean_posloaded_nrns, mean_negloaded_nrns, comp, dLP_idx, iLP_idx)
