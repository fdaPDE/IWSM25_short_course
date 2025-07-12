## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)


## ----echo=FALSE, out.width = '70%'--------------------------------------------
knitr::include_graphics("../data/fpca/HER2/hist_image.png")


## ----libraries, include=FALSE-------------------------------------------------
# Import the fdaPDE library
library(fdaPDE2)           # v. 2.0 (2025)
rm(list = ls())

# Load additional libraries and helper functions for plotting
source("../utils/graphics.R")


## ----warning=FALSE, include=FALSE---------------------------------------------
# Common plot settings
standard_plot_settings_fields <- function() {
  standard_plot_settings_fiends <- theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(angle = 45, hjust = 1)
    )
}


## ----data, fig.width=8, fig.height=5------------------------------------------
## [DATA]
path_data <- "../data/fpca/HER2/"
counts <- as.matrix(read.csv(paste(path_data,"counts.csv",sep=""),
                             row.names=1))
locations <- read.csv(file=paste(path_data,"locs.csv",sep=""), 
                      row.names=1)
# data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations
gene_names <- names(counts[,1])


## ----ERB22-gene---------------------------------------------------------------
idx.HER2 <- 22

plot_HER2 <- ggplot() +
    geom_tile(aes(x = locations[,1], y = locations[,2], 
                  fill = counts[idx.HER2,])) + 
    scale_fill_viridis_c(option = "magma") +
    coord_fixed() +
    standard_plot_settings_fields() +
    ggtitle(gene_names[idx.HER2])


## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("../data/fpca/HER2/reduced_hist_img.png")


## ----echo=FALSE---------------------------------------------------------------
plot_HER2


## ----mesh, results='hide'-----------------------------------------------------
## [MESH]
p <- pslg(
  P = locations
)
mesh_data <- triangulate(p, Y = FALSE, D = TRUE)
if (is.null(mesh_data$H)) mesh_data$H <- matrix(numeric(0), ncol = 2)
mesh_data <- triangulate(mesh_data, a = 0.1, q = 30, D = TRUE) ## ruppert's refinement
tissue <- triangulation(
  nodes    = mesh_data$P,
  cells    = mesh_data$T,
  boundary = mesh_data$PB
)


## ----mesh_plot, echo=FALSE, fig.width=6, fig.height=6-------------------------
plot(tissue)


## -----------------------------------------------------------------------------
## 1) Generate a rectangular grid
grid_step <- 1/6
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))

bbox <- SpatialPoints(locations)@bbox

xmin <- bbox[1,1]
xmax <- bbox[1,2]
ymin <- bbox[2,1]
ymax <- bbox[2,2]

x <- seq(xmin, xmax, by = grid_step)
y <- seq(ymin, ymax, by = grid_step)
grid <- expand.grid(x = x, y = y)


## 2) Adapt the rectangular grid to the convex hull of the locations
locations_sf <- st_as_sf(as.data.frame(locations), coords = c("x", "y"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = 4326)

# Compute convex hull
convex_hull <- st_convex_hull(st_union(st_as_sf(as.data.frame(locations), 
                                                coords = c("x", "y"), crs = 4326)))
# Keep only grid points within the convex hull
grid <- as.data.frame(st_coordinates(grid_sf[st_within(grid_sf, convex_hull, 
                                                               sparse = FALSE), ]))


## -----------------------------------------------------------------------------
# Firstly, construct a geoframe containing the sample means at locations
gf_mean <- geoframe(domain = tissue)
gf_mean$insert(layer = "sample_mean", 
          type = "point", 
          geo = locations, 
          data = data.frame(sample_mean=colMeans(counts)))

# Now let's smooth the colwise mean to obtain a smooth mean
f <- fe_function(tissue, type = "P1")
smooth_mean_m <- sr(formula = sample_mean ~ f, data = gf_mean)

sm_fit_log <- smooth_mean_m$fit(
  calibrator = gcv(
    optimizer = grid_search(grid = 10^seq(from=-5,to=-2,length.out=10)),
    edf = "hutch"
  )
)
# center the data
centred_counts <- sweep(counts, 2,
                        smooth_mean_m$fitted, 
                        FUN = "-")


## ----mean_visual,fig.width=6, fig.height=5------------------------------------
mean_function <- fe_function(domain = tissue, 
                 type = "P1", 
                 coeff = smooth_mean_m$f)
                     
ggplot() + 
      geom_tile(aes(x = grid[,1], y = grid[,2], fill = mean_function$eval(grid))) + 
      scale_fill_viridis_c(option = "magma") +
      coord_fixed() +
      standard_plot_settings_fields() +
      ggtitle("Smooth mean")


## -----------------------------------------------------------------------------
## construct a new geoframe
gf_fpca <- geoframe(domain = tissue)
gf_fpca$insert("gene_expression", type = "point", geo = locations)
gf_fpca[["gene_expression"]]$X <- t(centred_counts)

## modeling
fpca_m <- fpca("X", data = gf_fpca)
fpca_m$fit(
  npc = 6,
  calibrator = gcv(
    optimizer = grid_search(grid = 10^seq(from = -2, to = 2, by = 0.5))
  )
)


## ----variance-explained, fig.width=6, fig.height=5----------------------------
plot(apply(fpca_m$scores,2,var),type="b",
      xlab = "fPC index", ylab = "Variance explained")


## -----------------------------------------------------------------------------
plot_fPCs <- list()

fPCs.grid <- apply(fpca_m$pcs, 2, 
                   function(loadings){
                     f <- fe_function(domain = tissue, type = "P1", coeff = loadings)
                     f$eval(grid)
                   })

for (i in 1:3) {
  #flip fPCs to have coherent signs in the visualization
  if (max(fPCs.grid[, i],na.rm=T) < -min(fPCs.grid[, i],na.rm=T)) {
    fPCs.grid[, i] <- -fPCs.grid[, i]
  }
}

# plot using the same scale
global_min <- min(fPCs.grid[, 1:3],na.rm=T)
global_max <- max(fPCs.grid[, 1:3],na.rm=T)

plot_fPCs <- lapply(
  1:3, function(fPC_idx)
    ggplot() + 
      geom_tile(aes(x = grid[,1], y = grid[,2], fill = fPCs.grid[,fPC_idx])) + 
      scale_fill_viridis_c(option = "magma", limits = c(global_min, global_max)) +
      coord_fixed() +
      standard_plot_settings_fields() +
      ggtitle(paste("fPC",fPC_idx,sep=""))
  )


## ----echo=FALSE, out.width = '40%'--------------------------------------------
knitr::include_graphics("../data/fpca/HER2/reduced_hist_img.png")


## ----include=FALSE------------------------------------------------------------
library(patchwork)


## ----fig.width=8, fig.height=3------------------------------------------------
#plot the fPCs with a shared legend using the patchwork package
(plot_fPCs[[1]] + plot_fPCs[[2]] + plot_fPCs[[3]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
erb22_reconstruction <- 
  sweep(fpca_m$scores[,1:3] %*% t(fPCs.grid[,1:3]), 2, 
        mean_function$eval(grid), FUN = "+")[idx.HER2,]
  
plot_HER2_reconstruction <- ggplot() +
    geom_tile(aes(x = grid[,1], y = grid[,2], 
                  fill = erb22_reconstruction)) +
    scale_fill_viridis_c(option = "magma") + 
    coord_fixed() +
    standard_plot_settings_fields() +
    ggtitle("ERBB2 reconstruction")


## ----HER2_rec_plot, echo=FALSE, fig.width=6, fig.height=5---------------------
plot_HER2_reconstruction


## ----fig.width=6, fig.height=5------------------------------------------------
fpca_m$scores[idx.HER2,]

labels <- paste0("fPC", 1:3)
barplot(fpca_m$scores[idx.HER2,1:3],
        names.arg = labels, main = "ERB22 scores",
        ylim=c(-10,20))
abline(h=0)


## ----fig.width=6, fig.height=5------------------------------------------------
gene_names[which.max(fpca_m$scores[,1])]

#boxplot of scores1
boxplot(fpca_m$scores[,1], horizontal = T,
        main = "Score1")
points(max(fpca_m$scores[,1]), 1, col = "red", pch = 19, cex = 1.5)
text(max(fpca_m$scores[,1]), 1.1, labels = "ERB22", col = "red", pos = 2)

