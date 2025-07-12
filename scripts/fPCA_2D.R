## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  root.dir = "$(pwd)/../"
)


## ----fdaPDE, eval=TRUE--------------------------------------------------------
# Import the fdaPDE library
library(fdaPDE2)           # v. 2.0 (2025)
rm(list = ls())

# Load additional libraries and helper functions for plotting
source("../utils/graphics.R")


## ----include=FALSE------------------------------------------------------------
par(
  cex = 1.2,              # text size ~12pt
  font.main = 2,          # bold plot titles
  cex.main = 1.4,         # title size ~14pt
  col.main = "black",     # title color
  mar = c(5, 4, 4, 2) + 0.1,  # margins (bottom, left, top, right)
  las = 1                 # horizontal axis labels
)


## ----echo=FALSE, out.width = '70%'--------------------------------------------
knitr::include_graphics("../data/fPCA_2D/fPCA_2D_histological_image.png")


## ----data---------------------------------------------------------------------
## [DATA]
# Load the data
counts <- as.matrix(read.csv("../data/fPCA_2D/fPCA_2D_counts.csv", row.names = 1))
locations <- read.csv(file = "../data/fPCA_2D/fPCA_2D_locations.csv", row.names = 1)

# Data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations

# Gene names
gene_names <- names(counts[,1])


## ----ERB22-gene, fig.show='hide'----------------------------------------------
# ERB22-gene index
idx.HER2 <- 22

n_col <- 100
palette_ <- magma(n_col)
vals <- counts[idx.HER2,]
col <- palette_[as.numeric(cut(vals, breaks = n_col))]

par(oma = c(0, 0, 0, 0))
par(mar = c(0.5, 0.5, 0.5, 3))
plot(locations[, 1], locations[, 2], xlab = "", ylab = "", pch = 15,
     col = col, asp = 1, cex = 2, frame.plot=F, xaxt="n", yaxt="n")
image.plot(legend.only = TRUE, zlim = range(vals, na.rm = TRUE), col = palette_, 
           legend.args = list(side = 4), legend.lab = "", legend.mar = 2.25)


## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("../data/fPCA_2D/fPCA_2D_reduced_histological_image.png")


## ----data-plot, echo=FALSE----------------------------------------------------
# Static plot
par(oma = c(0, 0, 0, 0))
par(mar = c(0.5, 0.5, 0.5, 3))
plot(locations[, 1], locations[, 2], xlab = "", ylab = "", pch = 15,
     col = col, asp = 1, cex = 2, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(legend.only = TRUE, zlim = range(vals, na.rm = TRUE),
                   col = palette_, legend.args = list(side = 4),
                   legend.lab = "", legend.mar = 2.25)


## ----mesh-rtriangle, results='hide'-------------------------------------------
## [MESH]
# Create a planar straight line graph object
p <- pslg(P = locations)

# Create a regular mesh of the spatial domain
mesh <- triangulate(p, Y = FALSE, D = TRUE)
if (is.null(mesh$H)) mesh$H <- matrix(numeric(0), ncol = 2)

# Create a regular mesh of the spatial domain
mesh <- triangulate(mesh, a = 0.1, q = 30, D = TRUE)


## ----mesh-fdapde--------------------------------------------------------------
# Set up the triangulation for fdaPDE
mesh_tissue <- triangulation(
  nodes    = mesh$P,
  cells    = mesh$T,
  boundary = mesh$PB
)


## ----geoframe-----------------------------------------------------------------
# Create the geoframe
tissue <- geoframe(domain = mesh_tissue)
tissue


## ----geoframe-layer-----------------------------------------------------------
# Add layer with data to the geoframe object
tissue$insert(layer = "sample_mean", type = "point", 
              geo = locations, data = data.frame(sample_mean = colMeans(counts)))
tissue


## ----smoothing----------------------------------------------------------------
## [ISOTROPIC SMOOTHING WITH OPTIMAL SMOOTHING PARAMETER]
# Set up the finite element function (order 1)
f <- fe_function(domain = mesh_tissue, type = "P1")

# Proposed values for the smoothing parameter
lambda_grid <- 10^seq(from = -5, to = -2, length.out = 10)

# Isotropic smoothing model
smoothing_model <- sr(formula = sample_mean ~ f, data = tissue)

# Isotropic smoothing fit
smoothing_fit <- smoothing_model$fit(
  calibrator = gcv(
    optimizer = grid_search(grid = lambda_grid)
  )
)


## ----mean_visual, fig.width=5-------------------------------------------------
# Static plot
par(mar = c(0.5, 0.5, 0.5, 3), oma = c(0, 0, 0, 0), mfrow = c(1,1))
plot(f, palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(legend.only = TRUE, zlim = range(vals, na.rm = TRUE),
                   col = palette_, legend.args = list(side = 4),
                   legend.lab = "", legend.mar = 2.25)


## ----center-data--------------------------------------------------------------
# Center the data
centered_counts <- sweep(counts, 2, smoothing_model$fitted, FUN = "-")

# Create the geoframe
tissue_fpca <- geoframe(domain = mesh_tissue)

# Add layer with data to the geoframe object
tissue_fpca$insert("gene_expression", type = "point", geo = locations)
tissue_fpca[["gene_expression"]]$X <- t(centered_counts)


## ----fpca---------------------------------------------------------------------
## [FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS]
# fPCA model
fpca_model <- fpca(column = "X", data = tissue_fpca)

# Proposed values for the smoothing parameter
lambda_grid <- 10^seq(from = -2, to = 2, by = 0.5)

# fPCA fit
fpca_fit <- fpca_model$fit(
  npc = 6,
  calibrator = gcv(
    optimizer = grid_search(grid = lambda_grid)
  )
)


## ----variance-explained, fig.width=8.3, fig.height=5--------------------------
# Static plot
plot(apply(X = fpca_model$scores, MARGIN = 2, FUN = var), type = "b",
     xlab = "fPC index", ylab = "Explained variance")


## ----pcs----------------------------------------------------------------------
# Select the optimal number of PCs
npc <- 3

# fPCs
fPCs <- fpca_model$pcs[,1:npc]
for (i in 1:npc) {
  # Flip fPCs to have coherent signs in the visualization
  if (max(fPCs[, i], na.rm = T) < -min(fPCs[, i], na.rm = T)) {
    fPCs[, i] <- -fPCs[, i]
  }
}


## ----extracted-fpcs, fig.width=8.3, fig.height=5------------------------------
# Static plot
par(mfrow=c(1,3), mar = c(3, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 0))
plot(fe_function(domain = mesh_tissue, coeff=fPCs[,1], type = "P1"), 
     palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(legend.only = TRUE, zlim = range(fPCs[,1], na.rm = TRUE), 
                   col = palette_, legend.args = list(side = 4), 
                   legend.lab = "", horizontal = TRUE, legend.mar = 2.25)

plot(fe_function(domain = mesh_tissue, coeff=fPCs[,2], type = "P1"), 
     palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(legend.only = TRUE, zlim = range(fPCs[,2], na.rm = TRUE), 
                   col = palette_, legend.args = list(side = 4), 
                   legend.lab = "", horizontal = TRUE, legend.mar = 2.25)

plot(fe_function(domain = mesh_tissue,coeff=fPCs[,3], type = "P1"), 
     palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(legend.only = TRUE, zlim = range(fPCs[,3], na.rm = TRUE), 
                   col = palette_, legend.args = list(side = 4),
                   legend.lab = "", horizontal = TRUE, legend.mar = 2.25)


## ----erbreconstructiion, fig.width=5------------------------------------------
# Select fPCs associated with the ERBB2 gene
erb22 <- (fpca_model$scores[,1:3] %*% t(fPCs[,1:3]))[idx.HER2,]

# Static plot
par(mar = c(0.5, 0.5, 0.5, 3), oma = c(0, 0, 0, 0), mfrow = c(1,1))
plot(fe_function(domain = mesh_tissue, coeff=erb22, type = "P1"), 
     palette = magma, frame.plot = FALSE, xaxt = "n", yaxt = "n")
fields::image.plot(legend.only = TRUE, zlim = range(erb22, na.rm = T), 
                   col = palette_, legend.args = list(side = 4), 
                   legend.lab = "", legend.mar = 2.25)


## ----erbb2-fpc1, fig.width=8.3, fig.height=5----------------------------------
# Scores
fpca_model$scores[idx.HER2,]

# Static plot
labels <- paste0("fPC", 1:3)
barplot(fpca_model$scores[idx.HER2,1:3], names.arg = labels,
        main = "ERB22 scores", ylim = c(-10, 20))
abline(h = 0, lty = 2, lwd = 2)


## ----erbb2-scores, fig.width=8.3, fig.height=5--------------------------------
# Gene names
gene_names[which.max(fpca_model$scores[,1])]

# Static plot
boxplot(fpca_model$scores[,1], horizontal = TRUE, main = "Score1")
points(max(fpca_model$scores[,1]), 1, col = "red3", pch = 19, cex = 1.5)
text(max(fpca_model$scores[,1]), 1.1, labels = "ERB22", col = "red3", pos = 2)

