## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)


## ----library fdaPDE, eval=TRUE------------------------------------------------
# Import the fdaPDE library
library(fdaPDE2)           # v. 2.0 (2025)
rm(list = ls())

# Load additional libraries and helper functions for plotting
source("utils/graphics.R")



## ----theme, echo=FALSE--------------------------------------------------------
theme_set(theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_line()
  )
)


## ----geometry, fig.width=8, fig.height=5--------------------------------------
# Import the boundary nodes and segments of the domain of interest
boundary_nodes    <- read.table(file = "data/QSRPDE_2D/QSRPDE_2D_boundary_nodes.txt",
                                header = TRUE)
boundary_segments <- read.table(file = "data/QSRPDE_2D/QSRPDE_2D_boundary_segments.txt")
 
# Define a planar straight-line graph modeling the Switzerland boundary
p <- pslg(
  P = boundary_nodes,
  S = boundary_segments  
)
mesh <- triangulate(p, Y = FALSE, D = TRUE)
if (is.null(mesh$H)) mesh$H <- matrix(numeric(0), ncol = 2)
mesh <- triangulate(mesh, a = 0.0045, q = 30, D = TRUE) ## ruppert's refinement



## ----geometry_sf, fig.width=8, fig.height=5-----------------------------------
# create sf mesh object
mesh_sfc = st_as_sfc(mesh)
mesh_sf = st_as_sf(mesh_sfc, crs = 4326)
boundary_sf = st_boundary(st_union(mesh_sf))



## ----triangulation, fig.width=8, fig.height=5---------------------------------
switzerland = triangulation(nodes = mesh$P, cells = mesh$T, boundary = mesh$PB)


## ----mesh, fig.width=8.25, fig.height=5---------------------------------------
mapview(switzerland, crs = 4326, map.type = "CartoDB.Positron", 
          col.regions = "transparent", lwd = 1.25, legend = FALSE, layer = "mesh")


## ----data, fig.width=8, fig.height=5------------------------------------------
data <- read.table(file = "data/QSRPDE_2D/QSRPDE_2D_data.txt", header=TRUE)

gf <- geoframe(domain = switzerland)
gf$insert(layer = "rainfall", type = "point", geo = c("lon", "lat"), data = data)
gf


## ----mapview_data, fig.width=8.25, fig.height=5-------------------------------
# # Interactive plot
# mapview(gf[["rainfall"]], crs = 4326, 
#         color_palettes = list("mako", "mako")) 



map1 <- mapview(gf[["rainfall"]], crs = 4326, color_palettes = list("mako"), 
                varnames = "rainfall", na.color = "transparent",
                layer.name = "rainfall") +
  mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(gf[["rainfall"]], crs = 4326, color_palettes = list("mako"), 
                varnames = "log.rainfall", na.color = "transparent",
                layer.name = "rainfall") +
  mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2)




## ----estimate, results="hide"-------------------------------------------------
# physics modeling
f <- fe_function(switzerland, type = "P1")
K <- matrix(
    c(1.229618413471819, 1.001009926596135, 1.001009926596135, 1.628164356689574),
    nrow = 2, ncol = 2, byrow = TRUE
)
# modeling
m <- qsr(log.rainfall ~ f, data = gf, level = 0.5, penalty = fe_elliptic(K = K))

# fit
lambda_grid = 10^seq(from = -7, to = -3, by = 0.2)
fit = m$fit( calibrator = gcv( optimizer = grid_search(grid = lambda_grid ), seed=425) )



## ----grid_GCV, fig.width=8.3, fig.height=5------------------------------------
# GCV indices
gcv = fit$values  

# Optimal value selected for the smoothing parameter
lambda_opt_grid = fit$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(x = log10(lambda_grid), y = gcv, type = "b",
     lwd = 2, xlab = expression(log[10](lambda)), ylab = "GCV score")
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")



## ----mapview_log50, fig.width=8.3, fig.height=5-------------------------------
# Interactive plot
map_log50 <- mapview(f, crs = 4326, col.regions = mako,
                na.color = "transparent", layer.name = "rainfall", exp_transform=TRUE) +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)
map_log50



## ----mapview_50, fig.width=8.3, fig.height=5----------------------------------
# Interactive plot
map50 <- mapview(exp(f), crs = 4326, col.regions = mako,
                na.color = "transparent", layer.name = "rainfall", exp_transform=TRUE) +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)
map50



## ----estimate90, results="hide"-----------------------------------------------
# modeling
m <- qsr(log.rainfall ~ f, data = gf, level = 0.9, penalty = fe_elliptic(K = K))

# fit
lambda_grid = 10^seq(from = -7, to = -4, by = 0.2)
fit = m$fit( calibrator = gcv( optimizer = grid_search(grid = lambda_grid ), seed=425) )



## ----grid_GCV90, fig.width=8.3, fig.height=5----------------------------------
# GCV indices
gcv = fit$values  

# Optimal value selected for the smoothing parameter
lambda_opt_grid = fit$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(x = log10(lambda_grid), y = gcv, type = "b",
     lwd = 2, xlab = expression(log[10](lambda)), ylab = "GCV")
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")



## ----mapview_90, fig.width=8.3, fig.height=5----------------------------------
# Interactive plot
map90 <- mapview(f, crs = 4326, col.regions = mako,
                na.color = "transparent", layer.name = "rainfall", exp_transform=TRUE) +
         mapview(boundary_sf, col.regions = "transparent", alpha.regions = 0,
                col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map90


