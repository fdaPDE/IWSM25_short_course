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

######### GEOMETRY

# Load physical domain
boundary_nodes <- read.table(
  file = "../data/QSRPDE_2D/QSRPDE_2D_boundary_nodes.txt",
  header = TRUE
)
# Create boundary segments (consecutive boundary nodes are connected)
boundary_segments <- cbind(1:(nrow(boundary_nodes)-1), 2:nrow(boundary_nodes))
boundary_segments <- rbind(boundary_segments, c(nrow(boundary_nodes), 1))

# plot
boundary_nodes_sf <- st_as_sf(
  x = rbind(boundary_nodes, boundary_nodes[1,]),
  coords = c("lon", "lat"), crs = 4326
)
boundary_polygon <- st_polygon(list(st_coordinates(boundary_nodes_sf)))
boundary_sfc <- st_sfc(geometry = boundary_polygon, crs = 4326)
boundary_sf <- st_sf(geometry = boundary_sfc)

# Interactive plot
mapview(
  boundary_sf,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
  legend = FALSE, layer.name = "domain"
)

# meshing (using RTriangle)
p <- pslg(P = boundary_nodes, S = boundary_segments)
mesh <- triangulate(p = p, Y = FALSE, D = TRUE)
if (is.null(mesh$H)) mesh$H <- matrix(numeric(0), ncol = 2)
# mesh refinement
mesh <- triangulate(mesh, a = 0.0045, q = 30, D = TRUE)

mesh = triangulation(
  nodes = mesh$P, 
  cells = mesh$T, 
  boundary = mesh$PB
)

# plot
mapview(
  boundary_sf,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
  legend = FALSE, layer = "domain") +
mapview(
  mesh, 
  crs = 4326, col.regions = "transparent", lwd = 1.25, legend = FALSE, 
  layer = "mesh"
)

######### DATA

# create a geoframe object
switzerland <- geoframe(domain = mesh)
data <- read.table(file = "../data/QSRPDE_2D/QSRPDE_2D_data.txt", header = TRUE)
head(data)
switzerland$insert(
  layer = "rainfall", type = "point", geo = c("lon", "lat"),  data = data
)
switzerland

names(switzerland[["rainfall"]])              # Variable names
ncol(switzerland[["rainfall"]])               # Number of variables
head(gf_locations(switzerland[["rainfall"]])) # (lon, lat) of measurements

map1 <- 
  mapview(
    switzerland[["rainfall"]], varnames = "rainfall", 
    crs = 4326, color_palettes = list("mako"), na.color = "transparent",
    layer.name = "rainfall") +
  mapview(
    boundary_sf, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map2 <- 
  mapview(
    switzerland[["rainfall"]], varnames = "log.rainfall",
    crs = 4326, color_palettes = list("mako"), na.color = "transparent",
    layer.name = "log.rainfall") +
  mapview(
    boundary_sf, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map1, map2)

######### PHYSICS 
# Lf = -div(K grad(f))

K <- matrix(
    c(1.229618413471819, 1.001009926596135, 1.001009926596135, 1.628164356689574),
    nrow = 2, ncol = 2, byrow = TRUE
)

######### MODELING
# physics informed quantile regression
f <- fe_function(mesh, type = "P1")
# smoothing parameter proposals
lambda_grid = 10^seq(from = -7, to = -3, by = 0.2)

m <- qsr(
  log.rainfall ~ f, 
  data = switzerland, 
  level = 0.5, ## median
  penalty = fe_elliptic(K = K)
)
fit = m$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid), seed = 425)
)

# inspect GCV curve
gcv <- fit$values  
# Optimal value selected for the smoothing parameter
lambda_opt_grid <- fit$optimum
lambda_opt_grid
# Plot of the GCV curve
par(family = "serif")
plot(
  x = log10(lambda_grid), y = gcv, 
  type = "b", lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV"
)
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend(
  "bottomright", lty = 2, lwd = 2, col = "royalblue",
  legend = TeX("$\\log_{10}(\\lambda_{grid})")
)

## median field plot
map_log50 <- 
  mapview(
    f, 
    crs = 4326, col.regions = mako, na.color = "transparent", 
    layer.name = "log.rainfall") +
  mapview(
    boundary_sf, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map_50 <- 
  mapview(
    exp(f), 
    crs = 4326, col.regions = mako, na.color = "transparent", 
    layer.name = "rainfall") +
  mapview(
    boundary_sf, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map_log50, map_50)

# 90% quantile estimation
f <- fe_function(mesh, type = "P1")
lambda_grid = 10^seq(from = -7, to = -4, by = 0.2)
m <- qsr(
  log.rainfall ~ f, 
  data = switzerland, 
  level = 0.9, ## 90% quantile
  penalty = fe_elliptic(K = K)
)
fit = m$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid), seed = 425)
)

# inspect GCV curve
gcv <- fit$values  
# Optimal value selected for the smoothing parameter
lambda_opt_grid <- fit$optimum
lambda_opt_grid
# Plot of the GCV curve
par(family = "serif")
plot(
  x = log10(lambda_grid), y = gcv, 
  type = "b", lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV"
)
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend(
  "bottomright", lty = 2, lwd = 2, col = "royalblue",
  legend = TeX("$\\log_{10}(\\lambda_{grid})")
)

map_log90 <- 
  mapview(
    f, 
    crs = 4326, col.regions = mako, na.color = "transparent", 
    layer.name = "log.rainfall") +
  mapview(
    boundary_sf, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map_90 <- 
  mapview(
    exp(f), 
    crs = 4326, col.regions = mako, na.color = "transparent", 
    layer.name = "rainfall") +
  mapview(
    boundary_sf, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map_log90, map_90)
