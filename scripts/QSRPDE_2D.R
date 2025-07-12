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


## ----theme, echo=FALSE--------------------------------------------------------
theme_set(theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_line()
  )
)


## ----domain-------------------------------------------------------------------
## [SPATIAL DOMAIN]
# Load the boundary nodes
boundary_nodes <- read.table(file = "../data/QSRPDE_2D/QSRPDE_2D_boundary_nodes.txt",
                             header = TRUE)

# Create boundary segments (consecutive boundary nodes are connected)
boundary_segments <- cbind(1:(nrow(boundary_nodes)-1), 2:nrow(boundary_nodes))

# Close the loop (connect the last node to the first)
boundary_segments <- rbind(boundary_segments, c(nrow(boundary_nodes), 1))


## ----mapview_boundary, fig.width=8.3, fig.height=5----------------------------
# Convert boundary into a sf object
boundary_nodes_sf <- st_as_sf(x = rbind(boundary_nodes, boundary_nodes[1,]),
                              coords = c("lon", "lat"), crs = 4326)
boundary_polygon <- st_polygon(list(st_coordinates(boundary_nodes_sf)))
boundary_sfc <- st_sfc(geometry = boundary_polygon, crs = 4326)
boundary_sf <- st_sf(geometry = boundary_sfc)

# Interactive plot
mapview(boundary_sf,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
        legend = FALSE, layer.name = "domain")


## ----geometry, fig.width=8.3, fig.height=5------------------------------------
# Define a planar straight-line graph object
p <- pslg(P = boundary_nodes, S = boundary_segments)

# Create a regular mesh of the spatial domain
mesh <- triangulate(p = p, Y = FALSE, D = TRUE)
if (is.null(mesh$H)) mesh$H <- matrix(numeric(0), ncol = 2)

# Refine the mesh of the spatial domain
mesh <- triangulate(mesh, a = 0.0045, q = 30, D = TRUE)


## ----triangulation, fig.width=8.3, fig.height=5-------------------------------
# Set up the triangulation for fdaPDE
mesh = triangulation(nodes = mesh$P, cells = mesh$T, boundary = mesh$PB)

# Nodes coordinates
head(mesh$nodes)

# Number of nodes
mesh$n_nodes

# Edges
head(mesh$edges)

# Number of edges
mesh$n_edges

# Triangles
head(mesh$cells)

# Number of triangles
mesh$n_cells

# Bounding box
mesh$bbox


## ----mesh, fig.width=8.3, fig.height=5----------------------------------------
# Interactive plot
mapview(boundary_sf,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
        legend = FALSE, layer = "domain") +
  mapview(mesh, crs = 4326,
          col.regions = "transparent", lwd = 1.25, legend = FALSE, layer = "mesh")


## ----data, fig.width=8.3, fig.height=5----------------------------------------
switzerland <- geoframe(domain = mesh)


## ----load data----------------------------------------------------------------
## [DATA]
# Load the data
data <- read.table(file = "../data/QSRPDE_2D/QSRPDE_2D_data.txt", header = TRUE)
head(data)

# Add layer with data to the geoframe object
switzerland$insert(layer = "rainfall", type = "point", geo = c("lon", "lat"), 
                   data = data)
switzerland

# Variable names
names(switzerland[["rainfall"]])

# Number of variables
ncol(switzerland[["rainfall"]])

# Locations of measurement stations (lon, lat)
head(gf_locations(switzerland[["rainfall"]]))


## ----mapview_data, fig.width=8.3, fig.height=5--------------------------------
map1 <- mapview(switzerland[["rainfall"]], varnames = "rainfall", crs = 4326,
                color_palettes = list("mako"), na.color = "transparent",
                layer.name = "rainfall") +
  mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(switzerland[["rainfall"]], crs = 4326, varnames = "log.rainfall",
                color_palettes = list("mako"), na.color = "transparent",
                layer.name = "log.rainfall") +
  mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2)


## -----------------------------------------------------------------------------
## [DIFFUSION TENSOR]
K <- matrix(
    c(1.229618413471819, 1.001009926596135, 1.001009926596135, 1.628164356689574),
    nrow = 2, ncol = 2, byrow = TRUE
)


## ----estimate_50, results="hide"----------------------------------------------
## [PHYSICS-INFORMED QUANTILE REGRESSION]
# Set up the finite element function (order 1)
f <- fe_function(mesh, type = "P1")

# Proposed value for the smoothing parameter
lambda_grid = 10^seq(from = -7, to = -3, by = 0.2)

# Physics informed smoothing model
m <- qsr(log.rainfall ~ f, data = switzerland, level = 0.5, 
         penalty = fe_elliptic(K = K)
)

# Physics informed smoothing fit with grid search for GCV minimization
fit = m$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid),
                   seed = 425)
)


## ----grid_GCV_50, fig.width=8.3, fig.height=5---------------------------------
# GCV indices
gcv <- fit$values  

# Optimal value selected for the smoothing parameter
lambda_opt_grid <- fit$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(x = log10(lambda_grid), y = gcv, type = "b",
     lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV")
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend("bottomright", lty = 2, lwd = 2, col = "royalblue",
       legend = TeX("$\\log_{10}(\\lambda_{grid})"))


## ----mapview_50, fig.width=8.3, fig.height=5----------------------------------
# Interactive plot
map_log50 <- mapview(f, crs = 4326, col.regions = mako,
                     na.color = "transparent", layer.name = "log.rainfall") +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                 col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map_50 <- mapview(exp(f), crs = 4326, col.regions = mako,
                 na.color = "transparent", layer.name = "rainfall") +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map_log50, map_50)


## ----estimate_90, results="hide"----------------------------------------------
## [PHYSICS-INFORMED QUANTILE REGRESSION]
# Set up the finite element function (order 1)
f <- fe_function(mesh, type = "P1")

# Proposed value for the smoothing parameter
lambda_grid = 10^seq(from = -7, to = -4, by = 0.2)

# Physics informed smoothing model
m <- qsr(log.rainfall ~ f, data = switzerland, level = 0.9, 
         penalty = fe_elliptic(K = K)
)

# Physics informed smoothing fit with grid search for GCV minimization
fit = m$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid),
                   seed = 425)
)


## ----grid_GCV_90, fig.width=8.3, fig.height=5---------------------------------
# GCV indices
gcv <- fit$values  

# Optimal value selected for the smoothing parameter
lambda_opt_grid <- fit$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(x = log10(lambda_grid), y = gcv, type = "b",
     lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV")
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend("bottomright", lty = 2, lwd = 2, col = "royalblue",
       legend = TeX("$\\log_{10}(\\lambda_{grid})"))


## ----mapview_90, fig.width=8.3, fig.height=5----------------------------------
# Interactive plot
map_log90 <- mapview(f, crs = 4326, col.regions = mako,
                     na.color = "transparent", layer.name = "log.rainfall") +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                 col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map_90 <- mapview(exp(f), crs = 4326, col.regions = mako,
                 na.color = "transparent", layer.name = "rainfall") +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map_log90, map_90)

