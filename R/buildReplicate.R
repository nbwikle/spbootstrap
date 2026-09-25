### buildReplicate.R
### Turn one compact bootstrap replicate into a full spatial object.

buildReplicate <- function(data, replicate, coord.names = NULL) {
  # Build one bootstrap replicate of `data`, in the same class as `data`,
  # from its compact form. Resampled observations carry their new locations.
  #
  # Args:
  #   data       : the original data passed to the sampler (sf, sp, terra
  #                SpatRaster, raster Raster*, or data.frame/matrix)
  #   replicate  : one compact replicate, either
  #                  - list(idx, xy), one element of spbbPoint(..., output =
  #                    "indices") or spbbPolygon(..., output = "indices"); or
  #                  - an integer vector, one column of the index matrix from
  #                    spbbGrid() (e.g. b$bootstraps$bl.5[, k])
  #   coord.names: character(2), x and y column names; required only when
  #                `data` is a data.frame/matrix
  #
  # Returns the replicate as an object of the same class as `data`.
  #
  # Example (cluster jobs): generate compact replicates once and save them,
  # then build and analyse one replicate per job:
  #   reps <- spbbPolygon(wells, iowa, n.boot = 500, block.l = 30000,
  #                       output = "indices")
  #   saveRDS(reps$bl.30000, "reps.rds")
  #   ## in job k:
  #   r <- buildReplicate(wells, readRDS("reps.rds")[[k]])

  if (is.list(replicate) && !is.data.frame(replicate)) {
    if (!all(c("idx", "xy") %in% names(replicate))) {
      stop("A list `replicate` must have elements `idx` and `xy`, as returned with output = \"indices\".",
           call. = FALSE)
    }
    checkPointData(data)
    idx <- replicate$idx
    xy  <- as.matrix(replicate$xy)
    if (length(idx) != nrow(xy) || ncol(xy) != 2) {
      stop("`replicate$xy` must be a two-column matrix with one row per element of `replicate$idx`.",
           call. = FALSE)
    }
    checkIndices(idx, nObs(data))
    return(movePoints(data, idx, xy, coord.names))
  }

  if (is.numeric(replicate)) {
    n <- nObs(data)
    if (length(replicate) != n) {
      stop(sprintf(
        "A grid replicate must have one index per cell of `data` (%d), not %d. Use spbbGrid() indices, not point indices.",
        n, length(replicate)
      ), call. = FALSE)
    }
    checkIndices(replicate, n)
    if (!inherits(data, "sf") && (is.data.frame(data) || is.matrix(data)) &&
        (is.null(coord.names) || !all(coord.names %in% colnames(data)))) {
      stop("`data` is a data.frame/matrix, so `coord.names` must give the x and y column names.",
           call. = FALSE)
    }
    return(gridReplicate(data, as.integer(replicate), coord.names))
  }

  stop("`replicate` must be list(idx, xy) from spbbPoint()/spbbPolygon() or a column of spbbGrid() indices.",
       call. = FALSE)
}


nObs <- function(data) {
  # Number of observations (rows, features or cells) in `data`.
  if (inherits(data, "SpatRaster")) return(terra::ncell(data))
  if (inherits(data, "Raster")) return(raster::ncell(data))
  if (inherits(data, "sfc")) return(length(data))
  if (inherits(data, "Spatial")) return(length(data))
  nrow(data)
}


checkIndices <- function(idx, n) {
  if (anyNA(idx) || any(idx < 1) || any(idx > n) || any(idx != round(idx))) {
    stop(sprintf(
      "Replicate indices must be whole numbers between 1 and %d (the number of observations in `data`); was the replicate made from this data?",
      n
    ), call. = FALSE)
  }
  invisible(TRUE)
}
