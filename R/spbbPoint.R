### spbbPoint.R
### Block bootstrap for point process data observed in a rectangular box.

inBox <- function(pts, corner, b.l){
  which(pts[,1] >= corner[1] & pts[,1] <= corner[1] + b.l &
          pts[,2] >= corner[2] & pts[,2] <= corner[2] + b.l)
}


spbbPoint <- function(
  data, coord.names = NULL, n.boot, block.l, box.x = NULL, box.y = NULL,
  shift = TRUE, output = c("objects", "indices")
){
  # Create block bootstrap samples of a point pattern observed in a rectangular
  # box. Square source blocks are placed uniformly at random within the box,
  # and their points are translated onto a regular tiling of the box.
  #
  # Args:
  #   data       : point observations, as an sf/sfc object (POINT geometry),
  #                sp SpatialPoints* object, or data.frame/matrix
  #   coord.names: character(2), x and y column names; required only when
  #                `data` is a data.frame/matrix
  #   n.boot     : number of bootstrap replicates
  #   block.l    : vector of block side lengths, in coordinate units (any
  #                positive real value); each must be no larger than the
  #                shorter side of the box
  #   box.x      : numeric(2), x limits of the observation window; if NULL
  #                (default), the x range of the data's bounding box is used
  #   box.y      : numeric(2), y limits of the observation window; if NULL
  #                (default), the y range of the data's bounding box is used
  #   shift      : if TRUE (default), the tiling is moved by a random offset
  #                (uniform in [0, block.l) in x and y) for each replicate, so
  #                tile edges do not fall in the same places every time; if
  #                FALSE, tiles always start at the box's lower-left corner
  #   output     : "objects" (default) returns each replicate as an object of
  #                the same class as `data`; "indices" returns each replicate
  #                in compact form, list(idx = rows of `data`, xy = new
  #                coordinates), which is far smaller (e.g. for saving to disk
  #                and running replicates as separate cluster jobs). Turn a
  #                compact replicate into a full object with buildReplicate().
  #
  # The bounding box is the smallest rectangle containing the observed points,
  # so it is usually slightly smaller than the true study area; supply `box.x`
  # and `box.y` whenever the study area is known.
  #
  # Returns a named list (one element per block length), each a length-n.boot
  # list of bootstrap samples. Each sample has the same class as `data`
  # (data.frame/matrix, sf, sfc, or sp) and contains copies of the resampled
  # observations with their coordinates moved to their new locations.
  # With output = "indices", each sample is instead list(idx, xy).

  output <- match.arg(output)
  checkPointData(data)
  checkLonLat(list(data = data))
  reps <- pointReplicates(data, coord.names, n.boot, block.l, box.x, box.y, shift)
  if (output == "indices") return(reps)
  lapply(reps, lapply, function(r) movePoints(data, r$idx, r$xy, coord.names))
}


pointReplicates <- function(
  data, coord.names = NULL, n.boot, block.l, box.x = NULL, box.y = NULL, shift = TRUE
){
  # Workhorse for spbbPoint(): returns, for each block length, a length-n.boot
  # list of replicates in compact form, list(idx = rows of `data`, xy = new
  # coordinates), so that callers can build full objects one at a time.

  # grab data, determine box size
  checkPointData(data)
  pts <- spCoords(data, coord.names, arg = "data")

  box <- resolveBox(pts, box.x, box.y)
  box.x <- box$x; box.y <- box$y
  min.x <- min(box.x); max.x <- max(box.x); l.x <- max.x - min.x
  min.y <- min(box.y); max.y <- max(box.y); l.y <- max.y - min.y

  if (missing(block.l) || length(block.l) == 0) {
    stop("`block.l` must be supplied.", call. = FALSE)
  }
  if (!is.numeric(block.l) || anyNA(block.l) || any(block.l <= 0) ||
      any(block.l > min(l.x, l.y))) {
    stop(sprintf(
      "`block.l` must be greater than 0 and at most %g (the shorter side of the box).",
      min(l.x, l.y)
    ), call. = FALSE)
  }

  n.out <- sum(pts[, 1] < min.x | pts[, 1] > max.x | pts[, 2] < min.y | pts[, 2] > max.y)
  if (n.out > 0) {
    warning(sprintf(
      "%d of %d points lie outside `box.x` x `box.y` and will never be resampled.",
      n.out, nrow(pts)
    ), call. = FALSE)
  }

  # save bootstrap results
  bootstraps <- vector("list", length(block.l))

  for (m in seq_along(block.l)) {
    # block length
    b.l <- block.l[m]

    # number of blocks (one more in each direction when the tiling is shifted)
    n.tiles.x <- ceiling(l.x / b.l) + shift
    n.tiles.y <- ceiling(l.y / b.l) + shift
    checkTileCount(n.tiles.x * n.tiles.y, b.l)
    n.blocks  <- n.tiles.x * n.tiles.y

    # save bootstrap samples
    boot.samples <- list()

    for (k in 1:n.boot) {

      # tiling offset for this replicate (tiles start at the box's corner if 0)
      off <- if (shift) runif(2, 0, b.l) else c(0, 0)

      # sample corners
      corners.k <- cbind(
        runif(n = n.blocks, min = min.x, max = max.x - b.l),
        runif(n = n.blocks, min = min.y, max = max.y - b.l)
      )

      pts.lst <- list()
      pts.idx.lst <- list()

      for (i in 1:n.tiles.x){
        for (j in 1:n.tiles.y){

          # current index
          curr.idx <- (i - 1) * n.tiles.y + j

          # determine which points are in the box
          new.pts <- inBox(pts = pts, corner = corners.k[curr.idx,], b.l = b.l)
          orig.pts <- pts[new.pts, , drop = FALSE]

          # standardize the points
          std.pts <- cbind(
            orig.pts[, 1] - corners.k[curr.idx,1], orig.pts[, 2] - corners.k[curr.idx,2]
          )

          # move points into new tiled grid (tiles start at the box's lower-left
          # corner, moved down and left by the offset)
          pts.lst[[curr.idx]] <- cbind(
            min.x - off[1] + (i - 1) * b.l + std.pts[,1],
            min.y - off[2] + (j - 1) * b.l + std.pts[,2]
          )
          pts.idx.lst[[curr.idx]] <- new.pts
        }
      }

      # convert list to data frame
      pts.full <- do.call(rbind, pts.lst)
      pts.idx <- unlist(pts.idx.lst)

      # remove any points outside the boundary
      keep.pts <- which(pts.full[,1] <= max(box.x) & pts.full[,2] <= max(box.y) &
                          pts.full[,1] >= min(box.x) & pts.full[,2] >= min(box.y))
      pts.k <- pts.full[keep.pts, , drop = FALSE]
      pts.idx.k <- pts.idx[keep.pts]

      # save results
      boot.samples[[k]] <- list(idx = pts.idx.k, xy = pts.k)
    }
    bootstraps[[m]] <- boot.samples
  }

  names(bootstraps) <- paste0("bl.", block.l)
  bootstraps
}


movePoints <- function(data, idx, new.coords, coord.names = NULL) {
  # Subset `data` to rows `idx` and replace their coordinates with `new.coords`
  # (a two-column matrix), keeping the class of `data`.

  if (inherits(data, c("sf", "sfc"))) {
    geom <- sf::st_geometry(sf::st_as_sf(
      as.data.frame(new.coords), coords = 1:2, crs = sf::st_crs(data)
    ))
    if (inherits(data, "sfc")) return(geom)
    # subsetting the plain attribute table is much faster than `[.sf`
    gcol <- attr(data, "sf_column")
    out  <- sf::st_drop_geometry(data)[idx, , drop = FALSE]
    out[[gcol]] <- geom
    sf::st_sf(out, sf_column_name = gcol)

  } else if (inherits(data, "Spatial")) {
    crs <- methods::slot(data, "proj4string")
    colnames(new.coords) <- colnames(sp::coordinates(data))[1:2]
    if (methods::.hasSlot(data, "data")) {
      sp::SpatialPointsDataFrame(
        new.coords, data@data[idx, , drop = FALSE], proj4string = crs
      )
    } else {
      sp::SpatialPoints(new.coords, proj4string = crs)
    }

  } else {
    out <- data[idx, , drop = FALSE]
    out[, coord.names] <- new.coords
    out
  }
}
