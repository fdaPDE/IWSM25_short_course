# Import the fdaPDE library
library(fdaPDE2)           # v. 2.0 (2025)
rm(list = ls())

# Load additional libraries and helper functions for plotting
source("../utils/graphics.R")
theme_set(
  theme_minimal() +
    theme(
      axis.title = element_text(),
      axis.text = element_text(),
      axis.ticks = element_line(),
      panel.grid = element_line()
    )
)

# R-base plotting parameters
par(
  cex = 1.2,                  # text size ~12pt
  font.main = 2,              # bold plot titles
  cex.main = 1.4,             # title size ~14pt
  col.main = "black",         # title color
  mar = c(5, 4, 4, 2) + 0.1,  # margins (bottom, left, top, right)
  las = 1                     # horizontal axis labels
)

######### GEOMETRY

# Load physical domain
locations <- read.csv(
  file = "../data/fPCA_2D/fPCA_2D_locations.csv", row.names = 1
)

# triangulate the convex-hull of the locations
p <- pslg(P = locations)
mesh <- triangulate(p, Y = FALSE, D = TRUE)
if (is.null(mesh$H)) mesh$H <- matrix(numeric(0), ncol = 2)
mesh <- triangulate(mesh, a = 0.1, q = 30, D = TRUE)

mesh_tissue <- triangulation(
  nodes    = mesh$P,
  cells    = mesh$T,
  boundary = mesh$PB
)
plot(mesh_tissue)

######### DATA
counts <- as.matrix(
  read.csv("../data/fPCA_2D/fPCA_2D_counts.csv", row.names = 1)
)

# Data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations

# Gene names
gene_names <- names(counts[,1])
idx.HER2 <- 22 # ERB22-gene index

# plot ERB22-gene expression
n_col <- 100
palette_ <- magma(n_col)
vals <- counts[idx.HER2,]
col <- palette_[as.numeric(cut(vals, breaks = n_col))]

par(oma = c(0, 0, 0, 0))
par(mar = c(0.5, 0.5, 0.5, 3))
plot(
  locations[, 1], locations[, 2], 
  xlab = "", ylab = "", pch = 15, col = col, asp = 1, cex = 2, 
  frame.plot=F, xaxt="n", yaxt="n"
)
image.plot(
  legend.only = TRUE, zlim = range(vals, na.rm = TRUE), col = palette_, 
  legend.args = list(side = 4), legend.lab = "", legend.mar = 2.25
)

# create geoframe object
tissue <- geoframe(domain = mesh_tissue)
tissue
# add point layer for sample mean
tissue$insert(
  layer = "sample_mean", 
  type = "point", 
  geo = locations, 
  data = data.frame(sample_mean = colMeans(counts))
)
tissue

# recover smooth mean field
f <- fe_function(domain = mesh_tissue, type = "P1")
lambda_grid <- 10^seq(from = -5, to = -2, length.out = 10) # lambda proposals

m <- sr(sample_mean ~ f, data = tissue)

# Isotropic smoothing fit
r <- m$fit(
  calibrator = gcv(
    optimizer = grid_search(grid = lambda_grid)
  )
)

# smooth mean plot
par(mar = c(0.5, 0.5, 0.5, 3), oma = c(0, 0, 0, 0), mfrow = c(1,1))
plot(f, palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(
  legend.only = TRUE, zlim = range(vals, na.rm = TRUE),
  col = palette_, legend.args = list(side = 4),
  legend.lab = "", legend.mar = 2.25
)

######### MODELING
centered_counts <- sweep(counts, 2, m$fitted, FUN = "-")

# set-up geoframe
tissue_fpca <- geoframe(domain = mesh_tissue)
tissue_fpca$insert("gene_expression", type = "point", geo = locations)
tissue_fpca[["gene_expression"]]$X <- t(centered_counts)

# functional principal component analysis
fpca_model <- fpca(column = "X", data = tissue_fpca)

lambda_grid <- 10^seq(from = -2, to = 2, by = 0.5) # lambda proposals
fpca_fit <- fpca_model$fit(
  npc = 6,
  calibrator = gcv(
    optimizer = grid_search(grid = lambda_grid)
  )
)

## postprocess

## explained variance
plot(
  apply(X = fpca_model$scores, MARGIN = 2, FUN = var), type = "b",
  xlab = "fPC index", ylab = "Explained variance"
)

# as there is an elbow for npc <- 3, we select 3 as the optimal PC number
npc <- 3

# fPCs
fPCs <- fpca_model$pcs[,1:npc]
for (i in 1:npc) {
  # Flip fPCs to have coherent signs in the visualization
  if (max(fPCs[, i], na.rm = T) < -min(fPCs[, i], na.rm = T)) {
    fPCs[, i] <- -fPCs[, i]
  }
}

# PC functions plot
par(mfrow=c(1,3), mar = c(3, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 0))
plot(
  fe_function(domain = mesh_tissue, coeff=fPCs[,1], type = "P1"), 
  palette = magma, 
  frame.plot = FALSE, 
  xaxt = "n", 
  yaxt = "n"
)
fields::image.plot(
  legend.only = TRUE, zlim = range(fPCs[,1], na.rm = TRUE), 
  col = palette_, legend.args = list(side = 4), 
  legend.lab = "", horizontal = TRUE, legend.mar = 2.25
)

plot(
  fe_function(domain = mesh_tissue, coeff=fPCs[,2], type = "P1"), 
  palette = magma, 
  frame.plot = FALSE, 
  xaxt = "n", 
  yaxt = "n"
)
fields::image.plot(
  legend.only = TRUE, zlim = range(fPCs[,2], na.rm = TRUE), 
  col = palette_, legend.args = list(side = 4), 
  legend.lab = "", horizontal = TRUE, legend.mar = 2.25
)

plot(
  fe_function(domain = mesh_tissue, coeff=fPCs[,3], type = "P1"), 
  palette = magma, 
  frame.plot = FALSE, 
  xaxt = "n", 
  yaxt = "n"
)
fields::image.plot(
  legend.only = TRUE, zlim = range(fPCs[,3], na.rm = TRUE), 
  col = palette_, legend.args = list(side = 4), 
  legend.lab = "", horizontal = TRUE, legend.mar = 2.25
)

# use estimated PCs to reconstruct gene expression

# Select fPCs associated with the ERBB2 gene
erb22 <- (fpca_model$scores[,1:3] %*% t(fPCs[,1:3]))[idx.HER2,]

par(mar = c(0.5, 0.5, 0.5, 3), oma = c(0, 0, 0, 0), mfrow = c(1,1))
plot(
  fe_function(domain = mesh_tissue, coeff=erb22, type = "P1"), 
  palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n"
)
fields::image.plot(
  legend.only = TRUE, zlim = range(erb22, na.rm = T), col = palette_, 
  legend.args = list(side = 4), legend.lab = "", legend.mar = 2.25
)

# PC contributions to the expression of gene ERBB2
fpca_model$scores[idx.HER2,] # plot the scores

# Static plot
labels <- paste0("fPC", 1:3)
barplot(
  fpca_model$scores[idx.HER2,1:3], names.arg = labels, main = "ERB22 scores", 
  ylim = c(-10, 20)
)
abline(h = 0, lty = 2, lwd = 2)

# Gene names
gene_names[which.max(fpca_model$scores[,1])]

# Static plot
boxplot(fpca_model$scores[,1], horizontal = TRUE, main = "Score1")
points(max(fpca_model$scores[,1]), 1, col = "red3", pch = 19, cex = 1.5)
text  (max(fpca_model$scores[,1]), 1.1, labels = "ERB22", col = "red3", pos = 2)

