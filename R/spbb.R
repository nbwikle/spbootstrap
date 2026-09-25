### spbb.R
### Spatial block bootstrap for irregular (sf) point data.

spbb <- function(data, n.boot, block.size, boundary = NULL, resamp.size = NULL,
                 buffer = NULL, regular = FALSE, n.x = 0, n.y = 0){
  # Generates spatial block bootstrap samples for a given sf object
  #   Superseded by spbbPolygon(), which is faster and draws every source
  #   block uniformly from all squares inside the boundary; kept so earlier
  #   analyses can be reproduced.
  #   (Not optimal for a regular grid; use spbbGrid instead.)
  # Note: 'data' = sf object

  checkLonLat(list(data = data))

  ### set-up sf grid for sampling

  # create boundary (via convex hull)
  if (is.null(boundary)){
      boundary <- data %>% st_combine() %>%
        st_convex_hull()
  }

  if (!is.null(buffer)){
    boundary <- boundary %>%
      st_buffer(dist = 0.5, joinStyle  = "MITRE", mitreLimit = 2)
  }

  # create initial grid
  grid <- st_make_grid(
    boundary,
    n = block.size,
    what = "polygons",
    square = TRUE
  )

  # set resample size
  if (is.null(resamp.size)){
    resamp.size <- length(grid)
  }

  # determine side lengths of polygons
  side.length <- c(
    (st_bbox(grid)[3] - st_bbox(grid)[1]) / block.size[1],
    (st_bbox(grid)[4] - st_bbox(grid)[2]) / block.size[2]
  )

  if (regular){
    x.length <- n.x / block.size[1]
    y.length <- n.y / block.size[2]
  }

  ### select overlapping blocks

  boxes <- list()
  points <- list()
  df.only <- list()

  for (k in 1:resamp.size){

    if (regular){
      shift.k <- grid + c(
        sample(c(-1,1),1) * sample(0:(x.length-1), 1),
        sample(c(-1,1),1) * sample(0:(y.length-1), 1)
      )
    } else {
      shift.k <- grid +
        runif(n = 2, min = -1, max = 1) * side.length
      st_crs(shift.k) <- st_crs(grid)
    }

    interior.k <- st_contains(boundary, shift.k)[[1]]

    n.int <- length(interior.k)
    samp.k <- sample(n.int, size = 1)

    boxes[[k]] <- shift.k[interior.k[samp.k],]
    intersect.pts <- st_contains(boxes[[k]], st_geometry(data))[[1]]
    df.only[[k]] <- data[intersect.pts,] %>% st_drop_geometry()
    points[[k]] <- data[intersect.pts,] %>% st_geometry()
  }

  ### create bootstrapped data
  bootstraps <- list()
  for (j in 1:n.boot){

    bt.j <- spbootstrap(boundary, grid, boxes, points, df = df.only)

    if (regular){
      # reorder to match raster output
      #   -> start in upper left corner, left-to-right (i.e., by row)
      xy <- st_coordinates(bt.j)
      bt.j <- bt.j[order(xy[,"Y"], rev(xy[,"X"]), decreasing = TRUE),]
    }

    bootstraps[[j]] <- bt.j
  }

  ### return bootstraps
  return(bootstraps)
}


spbootstrap <- function(boundary, grid, boxes, points, df){

  # number of available blocks for bootstrap
  n.samps <- length(boxes)

  # polygons within boundary
  interior <- st_contains(boundary, grid)[[1]]
  # polygons along edge
  edge <- st_overlaps(boundary, grid)[[1]]

  # bootstrap interior data:
  bootstrap.int <- list()
  data.int <- list()

  for (j in 1:length(interior)){
    # sample one of the selected blocks
    boot.j <- sample(n.samps, 1)
    # determine shift to interior point
    bbox.j <- st_bbox(grid[interior[j],])
    bbox.boot <- st_bbox(boxes[[boot.j]])
    shift.j <- bbox.j[1:2] - bbox.boot[1:2]
    # save bootstrapped points
    bootstrap.int[[j]] <- points[[boot.j]] + shift.j
    st_crs(bootstrap.int[[j]]) <- st_crs(points[[boot.j]])

    # save bootstrapped data
    data.int[[j]] <- df[[boot.j]]
  }


  if (length(edge) > 0){
    # bootstrap edge data:
    bootstrap.edge <- list()
    data.edge <- list()

    for (j in 1:length(edge)){
      # sample one of the selected blocks
      boot.j <- sample(n.samps, 1)
      # determine shift to interior point
      bbox.j <- st_bbox(grid[edge[j],])
      bbox.boot <- st_bbox(boxes[[boot.j]])
      shift.j <- bbox.j[1:2] - bbox.boot[1:2]
      # save bootstrapped points
      bootstrap.edge[[j]] <- points[[boot.j]] + shift.j
      st_crs(bootstrap.edge[[j]]) <- st_crs(points[[boot.j]])

      data.edge[[j]] <- df[[boot.j]]
    }
  }

  # combine points into single sfc object
  interior.pts <- do.call(c, bootstrap.int)
  interior.df <- do.call(rbind, data.int)

  if (length(edge) > 0){
    exterior.pts <- do.call(c, bootstrap.edge)
    exterior.df <- do.call(rbind, data.edge)
    all.pts <- c(interior.pts, exterior.pts)
    all.df <- rbind(interior.df, exterior.df)
    all.sf <- st_as_sf(cbind(all.df, all.pts))
  } else {
    all.sf <- st_as_sf(cbind(interior.df, interior.pts))
  }

  # restrict to original spatial domain
  st_crs(all.sf) <- st_crs(boundary)
  original.domain <- st_intersects(boundary, all.sf)[[1]]
  bootstrap.sf <- all.sf[original.domain,]

  # return bootstrapped samples
  return(bootstrap.sf)
}


# ### Example
#
# ## make sf object
# ecuador.sf <- st_as_sf(ecuador, coords = c("x","y"))
#
# ## generate bootstrap samples
# bootstraps <- spbb(
#   data = ecuador.sf,
#   n.boot = 5,
#   block.size = c(3, 5),
#   boundary = NULL,
#   resamp.size = 100
# )
#
# ## compare bootstrap samples
# par(mfrow = c(2,2), mar = c(1,1,1,1))
# plot(st_geometry(bootstraps[[1]]), pch = 20, cex = 0.25)
# plot(st_geometry(bootstraps[[2]]), pch = 20, cex = 0.25)
# plot(st_geometry(bootstraps[[3]]), pch = 20, cex = 0.25)
# plot(st_geometry(bootstraps[[4]]), pch = 20, cex = 0.25)
#

