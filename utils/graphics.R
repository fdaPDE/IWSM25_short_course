# ADDITIONAL LIBRARIES
library(dplyr)
library(ggmap)
library(ggplot2)
library(latex2exp)
library(leafsync)
library(mapview)
library(raster)
library(RTriangle)
library(sf)
library(viridis)
library(fields)

# "ggforce", "htmlwidgets",
# "patchwork", "plotrix",
# "RColorBrewer", "simcausal", "sp"

# SET SEED
set.seed(23)

# HELPER FUNCTION: CONVERT A mesh.2D INTO A sfc --------------------------------
st_as_sfc.triangulation <- function(x, crs = NULL, ...){
  
  polygon_list = apply(x$T, MARGIN = 1, FUN = function(elem){
    st_cast(st_linestring(x$P[elem,]), to = "POLYGON")
  })
  
  if(is.null(crs)){crs = NA_crs_}
  
  mesh_sf = st_sfc(polygon_list, crs = crs)
  
  return(mesh_sf)
}

# HELPER FUNCTION: ORDER EDGES
order_edges <- function(edges) {
  ordered <- edges[1, ]
  remaining <- edges[-1, , drop = FALSE]
  
  while (nrow(remaining) > 0) {
    last <- ordered[length(ordered)]
    idx <- which(remaining[,1] == last | remaining[,2] == last)[1]
    
    if (is.na(idx)) break 
    next_edge <- remaining[idx, ]
    next_node <- if (next_edge[1] == last) next_edge[2] else next_edge[1]
    
    ordered <- c(ordered, next_node)
    remaining <- remaining[-idx, , drop = FALSE]
  }
  return(ordered)
}

# HELPER FUNCTION: CONVERT A mesh.2D INTO A sfc --------------------------------
st_as_sfc.triangulation_2_2 <- function(x, crs = NULL, ...){
  
  polygon_list = apply(x$cells, MARGIN = 1, FUN = function(elem){
    st_cast(st_linestring(x$nodes[elem,]), to = "POLYGON")
  })

  if(is.null(crs)){crs = NA_crs_}

  mesh_sf = st_sfc(polygon_list, crs = crs)

  return(mesh_sf)
}

# HELPER FUNCTION: CONVERT A geoframe INTO A sfc -------------------------------
st_as_sfc.gf_areal <- function(x, crs = NULL, ...){
  if(is.null(crs)){crs = NA_crs_}
  
  polygons <- gf_polygons(x)
  polygons_list <- vector("list", length(polygons))
  for(p in 1:length(polygons_list)){
    current <- polygons[[p]]
    current_edges <- order_edges(current$edges)
    current_edges <- c(current_edges, current_edges[1])
    current_coords <- current$nodes[current_edges,]
    polygons_list[[p]] <- st_polygon(list(current_coords))
  }
  polygons_sfc <- st_sfc(polygons_list, crs = crs)

  return(polygons_sfc)
}

# HELPER FUNCTION: COMPUTE AREAL EVALUATION ------------------------------------
eval_areal <- function(x, layer, crs = NULL){
  
  incidence_matrix = layer[["gf__ptr__"]]$incidence_matrix
  idx = which(incidence_matrix == 1, arr.ind = TRUE)
  markers = idx[order(idx[, "col"]), "row"]

  n_regions = nrow(layer)
  val = vector(mode = "numeric", length = n_regions)
  for(i in seq_len(n_regions)){
    incidence = function(cell) {
      cell_id = get_private(cell)$id_ + 1
      if(markers[cell_id] == i) return(TRUE)
      else return(FALSE)
    }
    get_private(layer[["gf__ptr__"]])$mesh_$mark_cells(marker = i, predicate = incidence)
    val[i] = x$integral(i)*100
  }
  
  if(is.null(crs)) crs = NA_crs_
  f = st_sf(data = val, geometry = st_as_sfc(layer), crs = crs)

  return(f)
}

# MAPVIEW WRAPPER --------------------------------------------------------------
mapview <- function(x, ..., res_lon = NULL, res_lat = NULL, varnames = NULL,
                    color_palettes = NULL, crs = NULL) {
  if (inherits(x, "triangulation")) {
    x_sfc = st_as_sfc(x, crs = crs)
    return(mapview::mapview(x_sfc, ...))
  }
  if (inherits(x, "triangulation_2_2")) {
    x_sfc = st_as_sfc(x, crs = crs)
    return(mapview::mapview(x_sfc, ...))
  }
  if (inherits(x, "gf_point")) {
    if(is.null(crs)) crs = NA_crs_
    if(is.null(varnames)) varnames = names(x)
    if(length(varnames) == 1){
      if(varnames == ""){
        df = gf_locations(x)
        df = as.data.frame(df)
        names(df) = c("lon", "lat")
        df$id = 1:nrow(df)
        df_sf = st_as_sf(x = df, coords = c("lon", "lat"), crs = crs)
        return(mapview::mapview(df_sf, legend = FALSE, ...))
      }
    }
    if(is.null(color_palettes)) color_palettes = list(mode = "character", length = length(varnames))
    if(length(color_palettes) == 1 & length(varnames) >= 1) {
      tmp = color_palettes[[1]]
      color_palettes = lapply(1:length(varnames), function(i) tmp)
    }
    varidxs = match(varnames, names(x))
    nvar = length(varnames)
    maps <- list()
    for(i in 1:nvar){
      j = varidxs[i]
      data = do.call("$", list(get("x"), names(x)[j]))
      data[is.infinite(data)] = NA
      df = data.frame(data = data,
                      lon = gf_locations(x)[,1],
                      lat = gf_locations(x)[,2])
      names(df)[1] = names(x)[j]
      df_sf = st_as_sf(x = df, coords = c("lon", "lat"), crs = crs)
      if(length(color_palettes[[i]]) == 1 & color_palettes[[i]][1] == "") color_palettes[[i]] = "viridis"
      maps[[i]] = mapview::mapview(df_sf, layer = names(x)[j],
                                   col.regions = ifelse(
                                     length(color_palettes[[i]]) == 1,
                                     match.fun(color_palettes[[i]]),
                                     colorRampPalette(color_palettes[[i]])),
                                   ...)
    }
    if(nvar == 1){
      return(maps[[1]])
    } else {
      ncol = ifelse(nvar %% 3, yes = 2, no = 3)
      return(sync(maps, ncol = ncol))
    }
  }
  if (inherits(x, "gf_areal")) {
    if(is.null(crs)) crs = NA_crs_
    if(is.null(varnames)) varnames = names(x)
    if(is.null(color_palettes)) color_palettes = vector(mode = "character", length = length(varnames))
    if(length(color_palettes) == 1 & length(varnames) >= 1) {
      tmp = color_palettes[[1]]
      color_palettes = lapply(1:length(varnames), function(i) tmp)
    }
    varidxs = match(varnames, names(x))
    nvar = length(varnames)
    maps <- list()
    for(i in 1:nvar){
      j = varidxs[i]
      data = do.call("$", list(get("x"), names(x)[j]))
      data[is.infinite(data)] = NA
      polygons_sfc <- st_as_sfc(x)
      df = as.data.frame(data)
      names(df)[1] = names(x)[j]
      df_sf = st_as_sf(x = df, geometry = polygons_sfc, crs = crs)
      if(length(color_palettes[[i]]) == 1 & color_palettes[[i]][1] == "") color_palettes[[i]] = "viridis"
      maps[[i]] = mapview::mapview(df_sf, layer = names(x)[j],
                                   col.regions = ifelse(
                                     length(color_palettes[[i]]) == 1,
                                     match.fun(color_palettes[[i]]),
                                     colorRampPalette(colors = color_palettes[[i]])), ...) 
    }
    if(nvar == 1){
      return(maps[[1]])
    } else {
      ncol = ifelse(nvar %% 3, yes = 2, no = 3)
      return(sync(maps, ncol = ncol))
    }
  }
  if (inherits(x, "plottable_function")) {
    if(is.null(res_lon)) res_lon = 250
    if(is.null(res_lat)) res_lat = 250
    
    bbox = x$geometry$bbox
    r = raster(nrows = res_lat, ncols = res_lon,
               xmn = bbox[1,1], xmx = bbox[2,1],
               ymn = bbox[1,2], ymx = bbox[2,2],
               crs = crs)
    grid = coordinates(r)
    values(r) = x$eval(grid)
    
    return(mapview::mapview(r, ...))
  }

  mapview::mapview(x, ..., crs = crs)
}

# fPCA
plot.plottable_function <- function(x, palette = NULL, ...) {
  n_col <- 100
  if (is.null(palette)) {
    palette_ <- colorRampPalette(colors = c("lightyellow", "darkred"))(n_col)
  } else if(is.function(palette)) {
    palette_ <- palette(n_col)
  } else { # user-defined palette
    palette_ <- palette_
  }
  
  nodes <- x$geometry$nodes
  x_grid <- seq(min(nodes[, 1]), max(nodes[, 1]), length.out = 250)
  y_grid <- seq(min(nodes[, 2]), max(nodes[, 2]), length.out = 250)
  xy_grid <- expand.grid(x_grid, y_grid)
  ## evaluate fe_function at fine grid
  vals <- x$eval(xy_grid)
  
  col <- palette_[as.numeric(cut(vals, breaks = n_col))]
  plot(
    xy_grid[, 1],
    xy_grid[, 2],
    xlab = "",
    ylab = "",
    pch = 15,
    col = col,
    asp = 1,
    cex = .6, ...
  )
}
