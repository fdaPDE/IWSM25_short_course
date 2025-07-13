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
domain <- st_read(
  dsn = "../data/SRPDE_2D/domain/SRPDE_2D_domain.shx", quiet = TRUE)
domain
mapview(
  domain,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
  legend = FALSE, layer.name = "domain"
)

# meshing
# define boundary nodes
boundary_nodes <- st_cast(x = domain, "POINT", crs = 4326)
boundary_nodes <- st_coordinates(x = boundary_nodes)
boundary_nodes <- data.frame(lon = boundary_nodes[,1], lat = boundary_nodes[,2])
# remove the last node (duplicate node)
boundary_nodes <- boundary_nodes[-nrow(boundary_nodes),]

# define boundary segments
boundary_segments <- cbind(1:nrow(boundary_nodes), c(2:nrow(boundary_nodes), 1))

# use RTriangle to generate a regular mesh
boundary_pslg <- pslg(P = boundary_nodes, S = boundary_segments)
mesh_regular <- triangulate(p = boundary_pslg)
if (is.null(mesh_regular$H)) {
  mesh_regular$H <- matrix(data = numeric(0), ncol = 2)
}
# mesh refinement
mesh_regular <- triangulate(p = mesh_regular, a = 0.01, q = 20, D = TRUE)

# wrap the triangulation with something understandable by fdaPDE
mesh <- triangulation(
  nodes = mesh_regular$P, 
  cells = mesh_regular$T,
  boundary = mesh_regular$PB
)

# you can query a set of informations
head(mesh$nodes) ## nodes coordinates
mesh$n_nodes     ## number of nodes
head(mesh$edges) ## edges' matrix
mesh$n_edges     ## number of edges
head(mesh$cells) ## triangles' matrix
mesh$n_cells     ## number of triangles
mesh$bbox        ## bounding box

mapview(
  domain,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
  legend = FALSE, layer = "domain") +
mapview(
  mesh, 
  crs = 4326, col.regions = "transparent", lwd = 1.25, legend = FALSE, 
  layer = "mesh"
)

######### DATA

# create a geoframe object
florida <- geoframe(domain = mesh)
florida

# load data as point layer
data <- read.table(file = "../data/SRPDE_2D/SRPDE_2D_data.txt")
florida$insert(
  layer = "ocean", type = "point", geo = c("lon", "lat"), data = data)
florida

names(florida[["ocean"]])              ## variable names
ncol(florida[["ocean"]])               ## number of variables
head(gf_locations(florida[["ocean"]])) ## measurement stations (lon, lat)

# Interactive plot
mapview(
  florida[["ocean"]], 
  crs = 4326,
  color_palettes = list(
    "inferno", "viridis", "viridis","viridis", "viridis", "viridis")
)

# load the SST data from NASA satellite images as a raster object
# we won't use these data during the fit, instead we use it for reference
load("../data/SRPDE_2D/raster/SST.RData")
SST
mapview(
  SST - 273.15, 
  col.regions = inferno, na.color = "transparent", layer.name = "SST [°C]") +
mapview(
  domain, 
  col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
  layer.name = "domain", legend = FALSE
)

######### MODELING
# let's try some physics-informed models for SST

## [ISOTROPIC SMOOTHING WITH FIXED SMOOTHING PARAMETER]
# set up the finite element function (order 1)
f_SST_iso_fixed <- fe_function(mesh, type = "P1")

lambda_fixed <- 0.00001 ## smoothing parameter proposal

## define the spatial regression model and fit
model_SST_iso_fixed <- sr(SST ~ f_SST_iso_fixed, data = florida)
fit_SST_iso_fixed   <- model_SST_iso_fixed$fit(lambda = lambda_fixed)

head(f_SST_iso_fixed$coeff)      # Fitted values at mesh nodes
head(model_SST_iso_fixed$fitted) # Fitted values at locations

# residuals at locations: response - fitted values
SST_iso_fixed_residuals <- 
  c(florida[["ocean"]]$SST - model_SST_iso_fixed$fitted)
summary(SST_iso_fixed_residuals)

map1 <- 
  mapview(
    f_SST_iso_fixed, ## fdaPDE integration with mapview
    crs = 4326, col.regions = inferno, na.color = "transparent", 
    layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map2 <- 
  mapview(
    SST - 273.15, 
    col.regions = inferno, na.color = "transparent", layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map1, map2)

# if you just want a rapid inspection of the results
plot(f_SST_iso_fixed, palette = inferno)

## [ISOTROPIC SMOOTHING WITH OPTIMAL SMOOTHING PARAMETER]
# set up the finite element function (order 1)
f_SST_iso_grid <- fe_function(mesh, type = "P1")

## smoothing parameter proposals
lambda_grid <- 10^seq(from = -6, to = -2, by = 0.2)

model_SST_iso_grid <- sr(formula = SST ~ f_SST_iso_grid, data = florida)

# fit with grid search for GCV minimization
fit_SST_iso_grid <- model_SST_iso_grid$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid))
)

head(f_SST_iso_grid$coeff)      ## Fitted values at mesh nodes
head(model_SST_iso_grid$fitted) ## Fitted values at locations

# residuals at locations: response - fitted values
SST_iso_grid_residuals <- c(florida[["ocean"]]$SST - model_SST_iso_grid$fitted)
summary(SST_iso_grid_residuals)

# optimal selected smoothing parameter
lambda_opt_grid <- fit_SST_iso_grid$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(
  x = log10(lambda_grid), y = fit_SST_iso_grid$values, 
  type = "b", lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV"
)
grid()
abline(v = log10(lambda_fixed),    lty = 2, lwd = 2, col = "lightblue3")
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue" )
legend("topleft", lty = c(2, 2), lwd = c(2, 2), 
       col = c("lightblue3", "royalblue"),
       legend = c(TeX("$\\log_{10}(\\lambda_{fixed})"),
                  TeX("$\\log_{10}(\\lambda_{grid})")))

# estimate plot
map1 <- 
  mapview(
    f_SST_iso_grid, 
    crs = 4326, col.regions = inferno, na.color = "transparent", 
    layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map2 <- 
  mapview(
    SST - 273.15, 
    col.regions = inferno, na.color = "transparent", layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map1, map2)

## isotropic smoothing doesn't seem a good modeling, let's try some more 
## fancy physics

## [TRANSPORT TERM]
# load the horizontal component of the transport term as a raster object
load("../data/SRPDE_2D/raster/Transport_x.RData")
Transport_x

# load the vertical component of the transport term as a raster object
load("../data/SRPDE_2D/raster/Transport_y.RData")
Transport_y

# get horizontal component of the transport term at (x,y)
get_transport_x <- function(x, y){
  value = extract(Transport_x, cbind(x, y))
  return(value)
}

# get vertical component of the transport term at (x,y)
get_transport_y <- function(x, y){
  value = extract(Transport_y, cbind(x, y))
  return(value)
}

# get the magnitude of the transport term at (x,y)
get_transport_coeff <- function(x, y){
  value = (get_transport_x(x = x, y = y))^2 + (get_transport_y(x = x, y = y))^2
  return(sqrt(value))
}

# plot transport field
df <- as.data.frame(coordinates(SST))
names(df) <- c("lon", "lat")
df$Transport_x <- ifelse(
  is.na(values(SST)), NA, get_transport_x(x = df$lon, y = df$lat))
df$Transport_y <- ifelse(
  is.na(values(SST)), NA, get_transport_y(x = df$lon, y = df$lat))
df$Transport_coeff <- ifelse(
  is.na(values(SST)), NA, get_transport_coeff(x = df$lon, y = df$lat))
# thin out the number of vectors
df <- df %>% filter(row_number() %% 1000 == 0)
# compute bounding box
bbox <- extent(SST)

## [STADIA MAPS API KEY]
# Insert here your API key; for link and instruction: ??register_stadiamaps
# register_stadiamaps(key = "--- your API key ---", write = TRUE)

# Map: run the three lines below once you have a working STADIA MAPS API key]
# map <- get_stadiamap(
#   bbox = c(
#     left = bbox@xmin, bottom = bbox@ymin, right = bbox@xmax, top = bbox@ymax),
#   zoom = 7, maptype = "alidade_smooth", color = "bw"
# )

# if you do not have a working STADIA MAPS API key, run the lines below
load(file = "../data/SRPDE_2D/raster/map.RData")

# Static plot
ggmap(map) +
  geom_segment(
    aes(xend = lon + Transport_x, yend = lat + Transport_y, 
        colour = Transport_coeff),
    data = df,
    arrow = arrow(angle = 23, length = unit(0.13, "inches")),
    size = 1.1,
    alpha = 1) +
  scale_color_gradientn(colours = rocket(n = 100, end = 0.87)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    text = element_text(size = 14),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)) +
  labs(color = "Velocity [m/s]") +
  guides(
    colour = guide_colourbar(
      barheight = unit(7, "cm"),barwidth = unit(0.55, "cm"), raster = TRUE, 
      ticks = FALSE)
    )

######### PHYSICS 
# Lf = -div(K grad(f)) + b * grad(f)

# diffusion tensor
K <- 0.0218 * diag(2)

# space-varying transport field
b <- function(points){
  output = array(data = 0, c(nrow(points),2))
  for (i in 1:nrow(points)){
    output[i,1] = get_transport_x(x = points[i,1], y = points[i,2])
    output[i,2] = get_transport_y(x = points[i,1], y = points[i,2])
  }
  return(output)
}

## [PHYSICS-INFORMED SMOOTHING WITH OPTIMAL SMOOTHING PARAMETER]
f_SST_physics <- fe_function(mesh, type = "P1")

# smoothing parameter proposals
lambda_grid <- 10^seq(from = -6, to = -2, by = 0.2)

# Physics-informed smoothing model
model_SST_physics <- sr(
  formula = SST ~ f_SST_physics, 
  data = florida,
  penalty = fe_elliptic(K = K, b = b)
)
fit_SST_physics <- model_SST_physics$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid))
)

# residuals at locations: response - fitted values
SST_physics_residuals <- c(florida[["ocean"]]$SST - model_SST_physics$fitted)
summary(SST_physics_residuals)

# residual boxplot
par(family = "serif")
boxplot(
  SST_iso_grid_residuals, SST_physics_residuals,
  col = c("forestgreen", "gold"), ylab = "residuals"
)
axis(
  1, at = 1:2, line = 0.5, tick = FALSE,
  labels = c("ISOTROPIC SR-PDE\n (GRID)", "PHYSICS-INFORMED SR-PDE\n (GRID)")
)

# optimal value selected for the smoothing parameter
lambda_opt_physics <- fit_SST_physics$optimum
lambda_opt_physics

map1 <- 
  mapview(
    florida[["ocean"]], varnames = "SST", 
    crs = 4326, layer.name = "SST [°C]", color_palettes = list("inferno"), 
    na.color = "transparent") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map2 <- 
  mapview(
    SST - 273.15, 
    col.regions = inferno, na.color = "transparent", layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map3 <- 
  mapview(
    f_SST_iso_grid, 
    crs = 4326, col.regions = inferno, na.color = "transparent", 
    layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map4 <- 
  mapview(
    f_SST_physics, 
    crs = 4326, col.regions = inferno, na.color = "transparent", 
    layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map1, map2, map3, map4, ncol = 2)

## modeling the quantity of dissolved oxigen

# Interactive plot
map1 <- 
  mapview(
    florida[["ocean"]], varnames = "SST", 
    crs = 4326, color_palettes = list("inferno"), na.color = "transparent",
    layer.name = "SST [°C]") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map2 <- 
  mapview(
    florida[["ocean"]], varnames = "DO", 
    crs = 4326, color_palettes = list("viridis"), na.color = "transparent",
    layer.name = "DO") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map1, map2)

# compute correlation
cor(x = florida[["ocean"]]$SST, y = florida[["ocean"]]$DO, method = "pearson")

######### PHYSICS 
# Lf = -div(K grad(f)) + b * grad(f)

# diffusion tensor
K <- 1.1341 * diag(2)

# space-varying transport field
b <- function(points){
  output = array(data = 0, c(nrow(points),2))
  for (i in 1:nrow(points)){
    output[i,1] = get_transport_x(x = points[i,1], y = points[i,2])
    output[i,2] = get_transport_y(x = points[i,1], y = points[i,2])
  }
  return(output)
}

## [ISOTROPIC SMOOTHING WITH A COVARIATE AND OPTIMAL SMOOTHING PARAMETER]
f_DO_physics <- fe_function(mesh, type = "P1")

model_DO_physics <- sr(
  DO ~ SST + f_DO_physics, 
  data = florida,
  penalty = fe_elliptic(K = K, b = b)
)
fit_DO_physics <- model_DO_physics$fit(
  calibrator = gcv(optimizer = grid_search(lambda_grid))
)

# residuals at locations: response - fitted values
DO_physics_residuals <- c(florida[["ocean"]]$DO - model_DO_physics$fitted)
summary(DO_physics_residuals)

# optimal value selected for the smoothing parameter
lambda_opt_physics <- fit_DO_physics$optimum
lambda_opt_physics

# estimate of beta parameter
model_DO_physics$beta

map1 <- 
  mapview(
    f_DO_physics, 
    crs = 4326, col.regions = viridis, na.color = "transparent", 
    layer.name = "DO") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

map2 <- 
  mapview(
    florida[["ocean"]], varnames = "DO", 
    crs = 4326, color_palettes = list("viridis"), na.color = "transparent",
    layer.name = "DO") +
  mapview(
    domain, 
    col.regions = "transparent", alpha.regions = 0, col = "black", lwd = 1.5, 
    layer.name = "domain", legend = FALSE
  )

sync(map1, map2)

