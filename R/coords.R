### coords.R
### Extract x/y coordinates from the spatial data types accepted by spbootstrap.

spCoords <- function(data, coord.names = NULL, arg = "data") {
  # Extract a two-column numeric data.frame of x/y coordinates from `data`.
  #
  # Args:
  #   data       : one of
  #                  - terra SpatRaster (all cells, in terra cell order)
  #                  - raster RasterLayer / RasterStack / RasterBrick
  #                  - sf / sfc object with POINT geometry
  #                  - sp Spatial* object (SpatialPoints*, SpatialPixels*,
  #                    SpatialGrid*)
  #                  - data.frame or matrix with coordinate columns
  #   coord.names: character(2), x and y column names; required for a
  #                data.frame/matrix, ignored for spatial classes
  #   arg        : argument name used in error messages
  #
  # Returns a data.frame with columns `x` and `y` (or `coord.names` for a
  # data.frame/matrix), one row per observation, in the same order as the
  # rows / cells / features of `data`.

  if (inherits(data, "SpatRaster")) {
    xy <- terra::xyFromCell(data, seq_len(terra::ncell(data)))

  } else if (inherits(data, "Raster")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required for raster objects.", call. = FALSE)
    }
    r  <- terra::rast(data)
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))

  } else if (inherits(data, c("sf", "sfc"))) {
    geom.type <- unique(as.character(sf::st_geometry_type(data)))
    if (!identical(geom.type, "POINT")) {
      stop(sprintf(
        "`%s` must have POINT geometry (found: %s); use sf::st_centroid() for polygon cells.",
        arg, paste(geom.type, collapse = ", ")
      ), call. = FALSE)
    }
    xy <- sf::st_coordinates(data)[, 1:2, drop = FALSE]

  } else if (inherits(data, "Spatial")) {
    if (!requireNamespace("sp", quietly = TRUE)) {
      stop("Package 'sp' is required for Spatial* objects.", call. = FALSE)
    }
    xy <- sp::coordinates(data)[, 1:2, drop = FALSE]

  } else if (is.data.frame(data) || is.matrix(data)) {
    if (is.null(coord.names) || length(coord.names) != 2) {
      stop(sprintf(
        "`%s` is a data.frame/matrix, so `coord.names` must give the x and y column names.",
        arg
      ), call. = FALSE)
    }
    missing.cols <- setdiff(coord.names, colnames(data))
    if (length(missing.cols) > 0) {
      stop(sprintf(
        "Column(s) not found in `%s`: %s", arg, paste(missing.cols, collapse = ", ")
      ), call. = FALSE)
    }
    return(as.data.frame(data)[, coord.names, drop = FALSE])

  } else {
    stop(sprintf(
      "`%s` must be a SpatRaster, Raster*, sf, Spatial*, data.frame or matrix (got %s).",
      arg, class(data)[1]
    ), call. = FALSE)
  }

  xy <- as.data.frame(xy)
  colnames(xy) <- c("x", "y")
  rownames(xy) <- NULL
  xy
}


resolveBox <- function(pts, box.x = NULL, box.y = NULL) {
  # Validate the observation window for point data, defaulting any missing
  # side to the bounding box of `pts` (same as sf::st_bbox / sp::bbox).
  #
  # Returns list(x = numeric(2), y = numeric(2)).

  default.box <- c(x = is.null(box.x), y = is.null(box.y))
  if (default.box["x"]) box.x <- range(pts[, 1])
  if (default.box["y"]) box.y <- range(pts[, 2])
  if (any(default.box)) {
    warning(sprintf(
      "No %s supplied; using the bounding box of the data (%s). Supplying the box is best when the study area is known.",
      paste0("`box.", names(default.box)[default.box], "`", collapse = " or "),
      paste(sprintf("%s: [%g, %g]", c("x", "y")[default.box],
                    c(box.x[1], box.y[1])[default.box],
                    c(box.x[2], box.y[2])[default.box]), collapse = ", ")
    ), call. = FALSE)
  }

  if (length(box.x) != 2 || length(box.y) != 2 || diff(range(box.x)) <= 0 ||
      diff(range(box.y)) <= 0) {
    stop("`box.x` and `box.y` must each be two distinct numbers.", call. = FALSE)
  }

  list(x = box.x, y = box.y)
}


isLonLat <- function(data) {
  # TRUE if `data` has a geographic (longitude/latitude) CRS; FALSE if it is
  # projected, has no CRS, or is a data.frame/matrix (whose CRS is unknown).
  out <- if (inherits(data, "SpatRaster")) {
    terra::is.lonlat(data, perhaps = FALSE, warn = FALSE)
  } else if (inherits(data, "Raster")) {
    raster::isLonLat(data)
  } else if (inherits(data, c("sf", "sfc"))) {
    sf::st_is_longlat(data)
  } else if (inherits(data, "Spatial")) {
    !sp::is.projected(data)
  } else {
    FALSE
  }
  isTRUE(out)
}


# Tracks whether an exported function further up the call stack has already
# checked its input, so nested calls (blockLength -> spboot -> spbbGrid) warn once.
.spbState <- new.env(parent = emptyenv())

checkLonLat <- function(objs, env = parent.frame()) {
  # Warn for each element of the named list `objs` that has unprojected
  # longitude/latitude coordinates. Only the outermost exported call checks;
  # the flag is cleared when that call returns.
  if (isTRUE(.spbState$checking)) return(invisible(FALSE))
  .spbState$checking <- TRUE
  do.call(on.exit, list(quote(.spbState$checking <- FALSE), add = TRUE), envir = env)

  lonlat <- vapply(objs, isLonLat, logical(1))
  for (arg in names(objs)[lonlat]) {
    warning(sprintf(paste(
      "`%s` has unprojected longitude/latitude coordinates, so block lengths and",
      "areas are in degrees and blocks are not square on the ground (a degree of",
      "longitude shrinks with latitude). Consider projecting to an equal-area or",
      "local coordinate system first, e.g. with sf::st_transform() or terra::project()."
    ), arg), call. = FALSE)
  }
  invisible(any(lonlat))
}


checkPointData <- function(data) {
  # Point-based samplers take sf/sp points or a data.frame/matrix, not rasters.
  if (inherits(data, c("SpatRaster", "Raster"))) {
    stop("`data` is a raster; use spbbGrid() for gridded data.", call. = FALSE)
  }
  if (inherits(data, "Spatial") && !inherits(data, "SpatialPoints")) {
    stop("sp input must be a SpatialPoints* object.", call. = FALSE)
  }
  invisible(TRUE)
}


checkTileCount <- function(n.tiles, b.l, max.tiles = 1e6) {
  # Guard against block lengths given in the wrong units (e.g. 1 when the
  # coordinates are in metres), which would create an enormous number of tiles.
  if (n.tiles > max.tiles) {
    stop(sprintf(paste(
      "`block.l` = %g would tile the study area with %.3g blocks. Block lengths are",
      "in the units of the coordinates (e.g. metres for UTM); is `block.l` in those units?"
    ), b.l, n.tiles), call. = FALSE)
  }
  invisible(TRUE)
}
