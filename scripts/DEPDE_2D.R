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
boundary_nodes <- read.table(file = "../data/DEPDE_2D/DEPDE_2D_boundary_nodes.txt")
boundary_segments <- read.table(file = "../data/DEPDE_2D/DEPDE_2D_boundary_segments.txt")

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
# Create a planar straight line graph object
boundary_pslg <- pslg(P = boundary_nodes, S = boundary_segments)
mesh <- triangulate(p = boundary_pslg)
if (is.null(mesh$H)) mesh$H <- matrix(data = numeric(0), ncol = 2)
mesh <- triangulate(p = mesh, a = 0.0001, q = 20, D = TRUE)

# create a triangulation object
mesh <- triangulation(
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
hampshire <- geoframe(domain = mesh)
hampshire
data <- read.table(file = "../data/DEPDE_2D/DEPDE_2D_data.txt")
head(data)
data <- data[1:2500,]

# add point pattern layer to the geoframe
hampshire$insert(
  layer = "diseases", type = "point", geo = c("lon", "lat"), data = data
)
hampshire

# plot
mapview(
  hampshire[["diseases"]], varnames = "", 
  crs = 4326, col.regions = "red3", cex = 2.5, layer.name = "locations") +
mapview(
  boundary_sf,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
  legend = FALSE, layer = "domain"
) 

######### MODELING
# point pattern estimation model

# smoothing parameter proposal
lambda_fixed <- 1e-2

model <- ppe(data = hampshire)
fit <- model$fit(
  lambda = lambda_fixed,
  optimizer = bfgs()
)

head(model$density)     # Density estimates at mesh nodes
head(model$log_density) # Log-density estimates at mesh nodes

# define a finite element function for the log_density field
f <- fe_function(domain =  mesh, type = "P1", coeff = model$log_density)

# plot
mapview(
  f, 
  crs = 4326, col.regions = inferno, na.color = "transparent", 
  layer.name = "ESTIMATE") +
mapview(
  boundary_sf,
  col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 2,
  legend = FALSE, layer.name = "domain"
)

