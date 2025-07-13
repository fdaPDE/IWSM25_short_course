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
data_sf <- st_read(dsn = "../data/GSRPDE_2D/GSRPDE_2D_data.shp", quiet = TRUE)
head(data_sf)
# number of data (i.e., of provinces)
n <- nrow(data_sf)

# extract provinces
provinces_sfc <- st_geometry(obj = data_sf)
provinces_sf <- st_as_sf(x = provinces_sfc)
# boundaries of each province
boundary_sf <- st_boundary(x = provinces_sf)

# Interactive plot
mapview(
  provinces_sf, 
  col.regions = "gray75", lwd = 1.5, legend = FALSE, layer.name = "provinces"
)

# meshing (imported from files)
mesh_nodes <- read.table(
  file = "../data/GSRPDE_2D/GSRPDE_2D_mesh_nodes.txt")
mesh_triangles <- read.table(
  file = "../data/GSRPDE_2D/GSRPDE_2D_mesh_triangles.txt")
mesh_nodesmarkers <- read.table(
  file = "../data/GSRPDE_2D/GSRPDE_2D_mesh_nodesmarkers.txt")

mesh <- triangulation(
  cells = mesh_triangles, 
  nodes = mesh_nodes,
  boundary = mesh_nodesmarkers
)

mapview(
  boundary_sf,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
  legend = FALSE, layer.name = "domain") +
  mapview(
    mesh, 
    crs = 4326, col.regions = "transparent", lwd = 1.25, legend = FALSE, 
    layer.name = "mesh"
  )

######### DATA

# create a geoframe object
italy <- geoframe(domain = mesh)
italy

# add layer with data to the geoframe object (directly from the shapefile)
italy$load_shp(
  layer = "fires", filename = "../data/GSRPDE_2D/GSRPDE_2D_data.shp"
)
italy

# you can query some informations from an areal layer
names(italy[["fires"]]) # Variable names
ncol(italy[["fires"]])  # Number of variables
head(gf_polygons(italy[["fires"]])[[1]]$nodes) # First province polygon nodes
head(gf_polygons(italy[["fires"]])[[1]]$nodes) # First province polygon edges

# Color palettes
color_palette_fire_counts <- c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026")
color_palette_population  <- c("#DEEBF7", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C")
color_palette_forest      <- c("#E5F5E0", "#A1D99B", "#74C476", "#31A354", "#006D2C")
color_palette_elevation   <- c("#006400", "#A2CD5A", "#F5DEB3", "#8B864E", "#8B2323")

mapview(
  italy[["fires"]], # fdaPDE types integration with mapview
  crs = 4326, varnames = c("FIRE_COUNT", "POP", "FOREST", "AVG_ELEV"),
  color_palettes = list(
    color_palette_fire_counts,
    color_palette_population,
    color_palette_forest,
    color_palette_elevation)
  )

######### MODELING

## [ISOTROPIC SMOOTHING WITH OPTIMAL SMOOTHING PARAMETER]
f_grid <- fe_function(mesh, type = "P1")

# smoothing parameter proposals
lambda_grid <- 10^seq(from = -6, to = 0, by = 0.25)

# poisson GLM for count data, isotropic smoothing
model_grid <- gsr(FIRE_COUNT ~ f_grid, data = italy, family = "poisson")
fit_grid <- model_grid$fit(
  calibrator = gcv(
    optimizer = grid_search(lambda_grid)
  )
)

head(f_grid$coeff)      # Fitted values at mesh nodes
head(model_grid$fitted) # Fitted values at locations
# residuals at locations: response - fitted values
grid_residuals <- c(italy[["fires"]]$FIRE_COUNT - model_grid$fitted)
summary(grid_residuals)

# inspect GCV curve
lambda_opt_grid <- fit_grid$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(
  x = log10(lambda_grid), y = fit_grid$values, 
  type = 'b', lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV"
)
grid()
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend(
  "topleft", lty = 2, lwd = 2, col = "royalblue", 
  legend = TeX("$\\log_{10}(\\lambda_{opt})")
)

# load locations of fires
locations <- st_read(
  dsn = "../data/GSRPDE_2D/GSRPDE_2D_locations.shx", quiet = TRUE)
map1 <- 
  mapview(
    f_grid, 
    crs = 4326, col.regions = color_palette_fire_counts, alpha.regions = 0.75, 
    na.color = "transparent", layer.name = "intensity") +
  mapview(
    boundary_sf, 
    color = "gray25", lwd = 1.5, legend = FALSE, layer.name = "provinces"
  )

map2 <- 
  mapview(
    locations, 
    legend = FALSE, col.regions = "black", alpha.regions = 0.75, stroke = FALSE, 
    cex = 2, layer.name = "locations") +
  mapview(
    boundary_sf, 
    color = "gray25", lwd = 1.5, legend = FALSE, layer.name = "provinces"
  )

sync(map1, map2)




# interactive plot
f_grid_areal <- eval_areal(x = f_grid, layer = italy[["fires"]], crs = 4326)
map1 <- 
  mapview(
    f_grid_areal, 
    col.regions = color_palette_fire_counts, na.color = "transparent", 
    layer.name = "ESTIMATE") +
   mapview(
     boundary_sf,
     col.regions = "transparent", alpha.regions = 0.25, col = "black", 
     lwd = 1.5, legend = FALSE, layer.name = "domain"
   )

map2 <- 
  mapview(
    italy[["fires"]], 
    crs = 4326, varnames = "FIRE_COUNT",
    color_palettes = list(color_palette_fire_counts)) +
  mapview(
    boundary_sf,
    col.regions = "transparent", alpha.regions = 0.25, col = "black", lwd = 1.5,
    legend = FALSE, layer.name = "domain"
  )

sync(map1, map2)
