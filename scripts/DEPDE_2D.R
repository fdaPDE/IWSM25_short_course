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


## ----spatial_domain-----------------------------------------------------------
## [SPATIAL DOMAIN]
# Load the boundary nodes
boundary_nodes <- read.table(file = "../data/DEPDE_2D/DEPDE_2D_boundary_nodes.txt")
head(boundary_nodes)

# Load the boundary segments
boundary_segments <- read.table(file = "../data/DEPDE_2D/DEPDE_2D_boundary_segments.txt")
head(boundary_segments)


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


## ----mesh---------------------------------------------------------------------
## [MESH]
# Create a planar straight line graph object
boundary_pslg <- pslg(P = boundary_nodes, S = boundary_segments)

# Create a regular mesh of the spatial domain
mesh <- triangulate(p = boundary_pslg)
if (is.null(mesh$H)) mesh$H <- matrix(data = numeric(0), ncol = 2)

# Refine the mesh of the spatial domain
mesh <- triangulate(p = mesh, a = 0.0001, q = 20, D = TRUE)
if (is.null(mesh$H)) mesh$H <- matrix(data = numeric(0), ncol = 2)

# Number of nodes
nrow(mesh$P)

# Number of triangles
nrow(mesh$T)


## ----triangulation------------------------------------------------------------
# Set up the triangulation for fdaPDE
mesh <- triangulation(nodes = mesh$P, cells = mesh$T,
                      boundary = mesh$PB)

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


## ----mapview_mesh, fig.width=8.3, fig.height=5--------------------------------
# Interactive plot
mapview(boundary_sf,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
        legend = FALSE, layer = "domain") +
  mapview(mesh, crs = 4326,
          col.regions = "transparent", lwd = 1.25, legend = FALSE, layer = "mesh")


## ----geoframe-----------------------------------------------------------------
# Create the geoframe
hampshire <- geoframe(domain = mesh)
hampshire


## ----data---------------------------------------------------------------------
## [DATA]
# Load the data
data <- read.table(file = "../data/DEPDE_2D/DEPDE_2D_data.txt")
head(data)

# Focus on cases occurred in 2001 (first 2500 cases)
data <- data[1:2500,]

# Add layer with data to the geoframe object
hampshire$insert(layer = "diseases", type = "point", geo = c("lon", "lat"), data = data)
hampshire


## ----mapview_data, fig.width=8.3, fig.height=5--------------------------------
# Interactive plot
mapview(hampshire[["diseases"]], varnames = "", crs = 4326, col.regions = "red3",
        cex = 2.5, layer.name = "locations") +
 mapview(boundary_sf,
         col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 1.5,
         legend = FALSE, layer = "domain") 


## ----density, results="hide"--------------------------------------------------
## [DENSITY ESTIMATION]
# Proposed value for the smoothing parameter
lambda_fixed <- 1e-2

# Density estimation model
model <- ppe(data = hampshire)

# Density estimation fit
fit <- model$fit(
  lambda = lambda_fixed,
  optimizer = bfgs()
)


## ----density_output-----------------------------------------------------------
# Density estimates at mesh nodes
head(model$density)

# Log-density estimates at mesh nodes
head(model$log_density)

# Finite element function
f <- fe_function(domain =  mesh, type = "P1", coeff = model$log_density)


## ----mapview_SST_iso_grid, fig.width=8.3, fig.height=5------------------------
# Interactive plot
mapview(f, crs = 4326, col.regions = inferno,
        na.color = "transparent", layer.name = "ESTIMATE") +
  mapview(boundary_sf,
        col.regions = "gray25", alpha.regions = 0.25, col = "black", lwd = 2,
        legend = FALSE, layer.name = "domain")

