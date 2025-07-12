## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)


## ----fdaPDE, eval=TRUE--------------------------------------------------------
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


## ----spatial_domain-----------------------------------------------------------
## [SPATIAL DOMAIN]
domain <- st_read(dsn = "data/SRPDE_2D/domain/SRPDE_2D_domain.shx")
domain


## ----mapview_spatial_domain, fig.width=8.3, fig.height=5----------------------
# Interactive plot
mapview(domain,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
        legend = FALSE, layer.name = "domain")


## ----mesh_regular-------------------------------------------------------------
## [MESH]
# Define boundary nodes
boundary_nodes <- st_cast(x = domain, "POINT", crs = 4326)
boundary_nodes <- st_coordinates(x = boundary_nodes)
boundary_nodes <- data.frame(lon = boundary_nodes[,1], lat = boundary_nodes[,2])

# Remove the last node (duplicate node)
boundary_nodes <- boundary_nodes[-nrow(boundary_nodes),]
head(boundary_nodes)

# Define boundary segments
boundary_segments <- cbind(1:nrow(boundary_nodes), c(2:nrow(boundary_nodes), 1))
head(boundary_segments)

# Create a planar straight line graph object
boundary_pslg <- pslg(P = boundary_nodes, S = boundary_segments)

# Create a regular mesh of the spatial domain
mesh_regular <- triangulate(p = boundary_pslg)
if (is.null(mesh_regular$H)) mesh_regular$H <- matrix(data = numeric(0), ncol = 2)

# Refine the mesh of the spatial domain
mesh_regular <- triangulate(p = mesh_regular, a = 0.01, q = 20, D = TRUE)
if (is.null(mesh_regular$H)) mesh_regular$H <- matrix(data = numeric(0), ncol = 2)

# Number of nodes
nrow(mesh_regular$P)

# Number of triangles
nrow(mesh_regular$T)


## ----triangulation------------------------------------------------------------
# Set up the triangulation for fdaPDE
mesh <- triangulation(nodes = mesh_regular$P, cells = mesh_regular$T,
                      boundary = mesh_regular$PB)

# Nodes coordinates
head(mesh$nodes)

# Number of nodes
mesh$n_nodes

# Edges
head(mesh$nodes)

# Number of edges
mesh$n_edges

# Triangles
head(mesh$cells)

# Number of triangles
mesh$n_cells

# Bounding box
mesh$bbox


## ----mapview_regular_mesh, fig.width=8.3, fig.height=5------------------------
# Interactive plot
mapview(domain,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
        legend = FALSE, layer = "domain") +
  mapview(mesh, crs = 4326,
          col.regions = "transparent", lwd = 1.25, legend = FALSE, layer = "mesh")


## ----geoframe-----------------------------------------------------------------
# Create the geoframe
florida <- geoframe(domain = mesh)
florida


## ----data---------------------------------------------------------------------
## [DATA]
# Load the data
data <- read.table(file = "data/SRPDE_2D/SRPDE_2D_data.txt")
head(data)

# Convert data into a sf object
data_sf <- st_as_sf(x = data, coords = c("lon", "lat"), crs = 4326)

# Add layer with data to the geoframe object
florida$insert(layer = "ocean", type = "point",
          geo = c("lon", "lat"), data = data)
florida


## ----mapview_data, fig.width=8.25, fig.height=5-------------------------------
# Interactive plot
mapview(florida[["ocean"]], crs = 4326,
        color_palettes = list("inferno", "viridis", "viridis", "viridis",
                              "viridis", "viridis"))


## ----SST----------------------------------------------------------------------
# Load the SST data from NASA satellite images as a raster object
load("data/SRPDE_2D/raster/SST.RData")
SST


## ----mapview_SST_raster, fig.width=8.3, fig.height=5--------------------------
# Interactive plot
mapview(SST - 273.15, col.regions = inferno, na.color = "transparent",
        layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)


## ----SST_iso_fixed, results="hide"--------------------------------------------
# Set up the finite element function
f_SST_iso_fixed <- fe_function(mesh, type = "P1")

# Proposed value for the smoothing parameter
lambda_fixed <- 0.00001

# Isotropic smoothing model
model_SST_iso_fixed <- sr(formula = SST ~ f_SST_iso_fixed, data = florida)

# Isotropic smoothing fit with fixed lambda
fit_SST_iso_fixed <- model_SST_iso_fixed$fit(lambda = lambda_fixed)


## ----SST_iso_fixed_output-----------------------------------------------------
# Fitted values at mesh nodes
head(f_SST_iso_fixed$coeff)

# Fitted values at locations
head(model_SST_iso_fixed$fitted)

# Residuals at locations: response - fitted values
SST_iso_fixed_residuals <- c(florida[["ocean"]]$SST - model_SST_iso_fixed$fitted)
summary(SST_iso_fixed_residuals)


## ----mapview_SST_iso_fixed, fig.width=8.3, fig.height=5-----------------------
# Interactive plot
map1 <- mapview(f_SST_iso_fixed, crs = 4326, col.regions = inferno,
                na.color = "transparent", layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(SST - 273.15, col.regions = inferno, na.color = "transparent",
                layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2)


## ----SST_iso_grid, results="hide"---------------------------------------------
# Proposed values for the smoothing parameter
lambda_grid <- 10^seq(from = -6, to = -2, by = 0.2)

# Set up the finite element function
f_SST_iso_grid <- fe_function(mesh, type = "P1")

# Isotropic smoothing model
model_SST_iso_grid <- sr(formula = SST ~ f_SST_iso_grid, data = florida)

# Isotropic smoothing fit with grid search for GCV minimization
fit_SST_iso_grid <- model_SST_iso_grid$fit(
  calibrator = gcv(optimizer = grid_search(grid = lambda_grid))
)


## ----SST_iso_grid_output------------------------------------------------------
# Fitted values at mesh nodes
head(f_SST_iso_grid$coeff)

# Fitted values at locations
head(model_SST_iso_grid$fitted)

# Residuals at locations: response - fitted values
SST_iso_grid_residuals <- c(florida[["ocean"]]$SST - model_SST_iso_grid$fitted)
summary(SST_iso_grid_residuals)


## ----SST_iso_grid_GCV, fig.width=8.3, fig.height=5----------------------------
# Optimal value selected for the smoothing parameter
lambda_opt_grid <- fit_SST_iso_grid$optimum
lambda_opt_grid

# Plot of the GCV curve
par(family = "serif")
plot(x = log10(lambda_grid), y = fit_SST_iso_grid$values, type = "b",
     lwd = 2, xlab = TeX("$\\log_{10}(\\lambda)$"), ylab = "GCV")
grid()
abline(v = log10(lambda_fixed), lty = 2, lwd = 2, col = "lightblue3")
abline(v = log10(lambda_opt_grid), lty = 2, lwd = 2, col = "royalblue")
legend("topleft", lty = c(2, 2), lwd = c(2, 2), col = c("lightblue3", "royalblue"),
       legend = c(TeX("$\\log(\\lambda_{fixed})"),
                  TeX("$\\log(\\lambda_{grid})")))


## ----mapview_SST_iso_grid, fig.width=8.3, fig.height=5------------------------
# Interactive plot
map1 <- mapview(f_SST_iso_grid, crs = 4326, col.regions = inferno,
                na.color = "transparent", layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(SST - 273.15, col.regions = inferno, na.color = "transparent",
                layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2)


## ----transport----------------------------------------------------------------
# Load the horizontal component of the transport term as a raster object
load("data/SRPDE_2D/raster/Transport_x.RData")
Transport_x

# Load the vertical component of the transport term as a raster object
load("data/SRPDE_2D/raster/Transport_y.RData")
Transport_y

# Get horizontal component of the transport term at (x,y)
get_transport_x <- function(x, y){
  value = extract(Transport_x, cbind(x, y))
  return(value)
}

# Get vertical component of the transport term at (x,y)
get_transport_y <- function(x, y){
  value = extract(Transport_y, cbind(x, y))
  return(value)
}

# Get the magnitude of the transport term at (x,y)
get_transport_coeff <- function(x, y){
  value = (get_transport_x(x = x, y = y))^2 + (get_transport_y(x = x, y = y))^2
  return(sqrt(value))
}


## ----include=FALSE------------------------------------------------------------
#api_key <- "..."
#register_stadiamaps(api_key)


## ----ggplot_transport, warning=FALSE, fig.width=8.3, fig.height=4-------------
# Create a data.frame with transport term values
df <- as.data.frame(coordinates(SST))
names(df) <- c("lon", "lat")
df$Transport_x <- ifelse(is.na(values(SST)), NA,
                         get_transport_x(x = df$lon, y = df$lat))
df$Transport_y <- ifelse(is.na(values(SST)), NA,
                         get_transport_y(x = df$lon, y = df$lat))
df$Transport_coeff <- ifelse(is.na(values(SST)), NA,
                             get_transport_coeff(x = df$lon, y = df$lat))

# Thin out the number of vectors
df <- df %>% filter(row_number() %% 1000 == 0)

## [STADIA MAPS API KEY]
# Insert here your API key; for link and instruction: ??register_stadiamaps
# register_stadiamaps(key = "---your API key---", write = TRUE)

# Compute the bounding box
bbox <- extent(SST)

# Static plot
ggmap(get_stadiamap(bbox = c(left = bbox@xmin, bottom = bbox@ymin,
                             right = bbox@xmax, top = bbox@ymax),
                    zoom = 7,
                    maptype = "alidade_smooth",
                    color = "bw")) +
  geom_segment(aes(xend = lon + Transport_x, yend = lat + Transport_y,
                   colour = Transport_coeff),
               data = df,
               arrow = arrow(angle = 23, length = unit(0.13, "inches")),
               size = 1.1,
               alpha = 1) +
  scale_color_gradientn(colours = rocket(n = 100, end = 0.87)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        text = element_text(size = 14)) +
  labs(color = "Velocity [m/s]") +
  theme(text = element_text(size = 14)) +
  guides(colour = guide_colourbar(barheight = unit(7, "cm"),
                                  barwidth = unit(0.55, "cm"),
                                  raster = TRUE, ticks = FALSE))


## ----PDE_parameters-----------------------------------------------------------
# Diffusion tensor
K <- 0.0218 * diag(2)

# Space-varying transport vector
b <- function(points){
  output = array(data = 0, c(nrow(points),2))
  for (i in 1:nrow(points)){
    output[i,1] = get_transport_x(x = points[i,1], y = points[i,2])
    output[i,2] = get_transport_y(x = points[i,1], y = points[i,2])
  }
  return(output)
}


## ----SST_physics, results="hide"----------------------------------------------
# Proposed values for the smoothing parameter
lambda_grid <- 10^seq(from = -6, to = -2, by = 0.2)

# Set up the finite element function
f_SST_physics <- fe_function(mesh, type = "P1")

# Physics-informed smoothing model
model_SST_physics <- sr(formula = SST ~ f_SST_physics, data = florida,
                        penalty = fe_elliptic(K = K, b = b))

# Physics-informed smoothing fit with Newton's method for GCV minimization
fit_SST_physics <- model_SST_physics$fit(calibrator = 
  gcv(optimizer = grid_search(grid = lambda_grid))
)


## ----SST_physics_output-------------------------------------------------------
# Fitted values at mesh nodes
head(f_SST_physics$coeff)

# Fitted values at locations
head(model_SST_physics$fitted)

# Residuals at locations: response - fitted values
SST_physics_residuals <- c(florida[["ocean"]]$SST - model_SST_physics$fitted)
summary(SST_physics_residuals)

# Optimal value selected for the smoothing parameter
lambda_opt_physics <- fit_SST_physics$optimum
lambda_opt_physics


## ----mapview_SST_physics, fig.width=8.3, fig.height=5-------------------------
# Interactive plot
map1 <- mapview(florida[["ocean"]], varnames = "SST", crs = 4326, layer.name = "SST [°C]",
                color_palettes = list("inferno"), na.color = "transparent") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(SST - 273.15, col.regions = inferno, na.color = "transparent",
                layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map3 <- mapview(f_SST_iso_grid, crs = 4326, col.regions = inferno,
                na.color = "transparent", layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map4 <- mapview(f_SST_physics, crs = 4326, col.regions = inferno,
                na.color = "transparent", layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2, map3, map4, ncol = 2)


## ----SST_boxplot, fig.width=8.3, fig.height=5---------------------------------
# Boxplot
par(family = "serif")
boxplot(SST_iso_grid_residuals, SST_physics_residuals,
        col = rocket(2), ylab = "residuals")
axis(1, at = 1:2, line = 0.5, tick = FALSE,
     labels = c("ISOTROPIC SR-PDE\n (GRID)", "PHYSICS-INFORMED SR-PDE\n (NEWTON)"))


## ----mapview_DO_data, fig.width=8.3, fig.height=5-----------------------------
# Interactive plot
map1 <- mapview(florida[["ocean"]], varnames = "SST", crs = 4326,
                color_palettes = list("inferno"), na.color = "transparent",
                layer.name = "SST [°C]") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(florida[["ocean"]], varnames = "DO", crs = 4326,
                color_palettes = list("viridis"), na.color = "transparent",
                layer.name = "DO") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2)


## ----correlation--------------------------------------------------------------
cor(x = florida[["ocean"]]$SST, y = florida[["ocean"]]$DO, method = "pearson")


## ----PDE_diffusion------------------------------------------------------------
# Diffusion tensor
K <- 1.1341 * diag(2)


## ----DO_physics---------------------------------------------------------------
# Set up the finite element function
f_DO_physics <- fe_function(mesh, type = "P1")

# Physics-informed smoothing model
model_DO_physics <- sr(formula = DO ~ SST + f_DO_physics, data = florida,
                        penalty = fe_elliptic(K = K, b = b))

# Physics-informed smoothing fit Grid Search method for GCV minimization
fit_DO_physics <- model_DO_physics$fit(
  calibrator = gcv(optimizer = grid_search(lambda_grid))
)


## ----DO_physics_output--------------------------------------------------------
# Fitted values at mesh nodes
head(f_DO_physics$coeff)

# Fitted values at locations
head(model_DO_physics$fitted)

# Residuals at locations: response - fitted values
DO_physics_residuals <- c(florida[["ocean"]]$DO - model_DO_physics$fitted)
summary(DO_physics_residuals)

# Optimal value selected for the smoothing parameter
lambda_opt_physics <- fit_DO_physics$optimum
lambda_opt_physics

# Estimate of parameter beta
model_DO_physics$beta


## ----mapview_DO_physics, fig.width=8.3, fig.height=5--------------------------
# Interactive plot
map1 <- mapview(f_DO_physics, crs = 4326, col.regions = viridis,
                na.color = "transparent", layer.name = "DO") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

map2 <- mapview(florida[["ocean"]], varnames = "DO", crs = 4326,
                color_palettes = list("viridis"), na.color = "transparent",
                layer.name = "DO") +
  mapview(domain, col.regions = "transparent", alpha.regions = 0,
          col = "black", lwd = 1.5, layer.name = "domain", legend = FALSE)

sync(map1, map2)

