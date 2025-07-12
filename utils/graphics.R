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
library(patchwork)
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
