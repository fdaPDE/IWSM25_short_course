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


## ----data---------------------------------------------------------------------
## [DATA]
# Load the data
data_sf <- st_read(dsn = "../data/GSRPDE_2D/GSRPDE_2D_data.shp", quiet = TRUE)
head(data_sf)

# Number of data (i.e., of provinces)
n <- nrow(data_sf)


## ----mapview_boundaries, fig.width=8.3, fig.height=5--------------------------
## [SPATIAL DOMAIN]
# Provinces
provinces_sfc <- st_geometry(obj = data_sf)
provinces_sf <- st_as_sf(x = provinces_sfc)
provinces_sf

# Boundaries of each province
boundary_sf <- st_boundary(x = provinces_sf)
boundary_sf

# Interactive plot
mapview(provinces_sf, col.regions = "gray75", lwd = 1.5,
        legend = FALSE, layer.name = "provinces")


## ----mesh---------------------------------------------------------------------
## [MESH]
# Load the mesh nodes
mesh_nodes <- read.table(file = "../data/GSRPDE_2D/GSRPDE_2D_mesh_nodes.txt")
head(mesh_nodes)

# Load the nodes markers
mesh_nodesmarkers <- read.table(file = "../data/GSRPDE_2D/GSRPDE_2D_mesh_nodesmarkers.txt")
head(mesh_nodesmarkers)

# Load the mesh triangles
mesh_triangles <- read.table(file = "../data/GSRPDE_2D/GSRPDE_2D_mesh_triangles.txt")
head(mesh_triangles)

# Set up the triangulation for fdaPDE
mesh <- triangulation(cells = mesh_triangles, nodes = mesh_nodes,
                      boundary = mesh_nodesmarkers)

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

# Bounding Box
mesh$bbox


## ----mapview_mesh, fig.width=8.3, fig.height=5--------------------------------
# Interactive plot
mapview(boundary_sf,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
        legend = FALSE, layer.name = "domain") +
  mapview(mesh, crs = 4326, col.regions = "transparent", lwd = 1.25,
          legend = FALSE, layer.name = "mesh")


## ----geoframe-----------------------------------------------------------------
# Create the geoframe
italy <- geoframe(domain = mesh)
italy


## ----geoframe_layer-----------------------------------------------------------
## [DATA]
# Add layer with data to the geoframe object (directly from the shapefile)
italy$load_shp(layer = "fires", filename = "../data/GSRPDE_2D/GSRPDE_2D_data.shp")
italy

# Variable names
names(italy[["fires"]])

# Number of variables
ncol(italy[["fires"]])

# First province polygon nodes
head(gf_polygons(italy[["fires"]])[[1]]$nodes)

# First province polygon edges
head(gf_polygons(italy[["fires"]])[[1]]$nodes)


## ----mapview_data, fig.width=8.3, fig.height=5--------------------------------
# Color palettes
color_palette_fire_counts <- c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026")
color_palette_population <- c("#DEEBF7", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C")
color_palette_forest <- c("#E5F5E0", "#A1D99B", "#74C476", "#31A354", "#006D2C")
color_palette_elevation <- c("#006400", "#A2CD5A", "#F5DEB3", "#8B864E", "#8B2323")

# Interactive plot
mapview(italy[["fires"]], crs = 4326, varnames = c("FIRE_COUNT", "POP", "FOREST", "AVG_ELEV"),
        color_palettes = list(color_palette_fire_counts,
                              color_palette_population,
                              color_palette_forest,
                              color_palette_elevation))


## ----solution_grid, results="hide"--------------------------------------------
## [ISOTROPIC SMOOTHING WITH OPTIMAL SMOOTHING PARAMETER]
# Set up the finite element function (order 1)
f_grid <- fe_function(mesh, type = "P1")

# Proposed values for the smoothing parameter
lambda_grid <- 10^seq(from = -6, to = 0, by = 0.25)

# Isotropic smoothing model
model_grid <- gsr(formula = FIRE_COUNT ~ f_grid, data = italy, family = "poisson")

# Isotropic smoothing fit with fixed lambda
fit_grid <- model_grid$fit(
  calibrator = gcv(
    optimizer = grid_search(lambda_grid)
  )
)


## ----grid_output--------------------------------------------------------------
# Fitted values at mesh nodes
head(f_grid$coeff)

# Fitted values at locations
head(model_grid$fitted)

# Residuals at locations: response - fitted values
grid_residuals <- c(italy[["fires"]]$FIRE_COUNT - model_grid$fitted)
summary(grid_residuals)


## ----solution_GCV, fig.width=8.3, fig.height=5--------------------------------
# Optimal value selected for the smoothing parameter
lambda_opt_grid <- fit_grid$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(x = log10(lambda_grid), y = fit_grid$values, type = 'b',
     lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV")
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend("topleft", lty = 2, lwd = 2, col = "royalblue",
       legend = TeX("$\\log_{10}(\\lambda_{opt})"))


## ----mapview_solution_grid_fit, fig.width=8.3, fig.height=5-------------------
# Compute areal estimate by province
f_grid_areal <- eval_areal(x = f_grid, layer = italy[["fires"]], crs = 4326)

# Interactive plot
map1 <- mapview(f_grid_areal, col.regions = color_palette_fire_counts,
                na.color = "transparent", layer.name = "ESTIMATE") +
 mapview(boundary_sf,
         col.regions = "transparent", alpha.regions = 0.25, col = "black", lwd = 1.5,
         legend = FALSE, layer.name = "domain")

map2 <- mapview(italy[["fires"]], crs = 4326, varnames = "FIRE_COUNT",
                color_palettes = list(color_palette_fire_counts)) +
 mapview(boundary_sf,
         col.regions = "transparent", alpha.regions = 0.25, col = "black", lwd = 1.5,
         legend = FALSE, layer.name = "domain")

sync(map1, map2)


## ----mapview_solution_point, fig.width=8.3, fig.height=5----------------------
# Load locations of fires
locations <- st_read(dsn = "../data/GSRPDE_2D/GSRPDE_2D_locations.shx", quiet = TRUE)

# Interactive plot
map1 <- mapview(f_grid, crs = 4326, col.regions = color_palette_fire_counts,
                alpha.regions = 0.75, na.color = "transparent",
                layer.name = "intensity") +
  mapview(boundary_sf, color = "gray25", lwd = 1.5,
          legend = FALSE, layer.name = "provinces")

map2 <- mapview(locations, legend = FALSE, col.regions = "black",
                alpha.regions = 0.75, stroke = FALSE, cex = 2,
                layer.name = "locations") +
  mapview(boundary_sf, color = "gray25", lwd = 1.5,
          legend = FALSE, layer.name = "provinces")

sync(map1, map2)

