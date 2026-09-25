### spsubGrid.R
### Window subsampling for gridded raster data.

spsubGrid <- function(
  data, coord.names, block.l, kappa = 0,
  pt.data = NULL, pt.coord.names = NULL
) {
  # Place overlapping subsampling windows (tiles) over a rectangular grid and
  # return the original-data indices of the points in each tile.
  # Optionally, also returns the indices of point process observations
  # (not on the grid) falling within each tile.
  #
  # Args:
  #   data          : data.frame with gridded observations
  #   coord.names   : character(2), column names for the x and y coordinates
  #   block.l       : integer (or integer vector), tile side length in grid cells
  #   kappa         : overlap proportion in [0, 1).
  #                     kappa = 0  -> non-overlapping tiles, step = block.l
  #                     kappa > 0  -> step = max(1, round(block.l * (1 - kappa)))
  #                   so neighbouring tiles share kappa * block.l cells on each side.
  #   pt.data       : (optional) data.frame of point process observations
  #   pt.coord.names: (optional) character(2), column names for x and y in pt.data
  #
  # Returns a list with:
  #   tiles  : named list (one element per block length) each containing
  #     indices    : (block.l^2) x n.tiles integer matrix of row indices into data,
  #                  one column per tile, in column-major (y-first) order within tile
  #     pt.indices : (if pt.data provided) length-n.tiles list of integer vectors,
  #                  each giving the row indices into pt.data for that tile
  #     n.tiles    : total number of tiles
  #     n.tiles.x  : number of tiles in the x direction
  #     n.tiles.y  : number of tiles in the y direction
  #     step       : actual step size used (grid cells)
  #     x.starts   : x-direction tile origin positions (standardized, 1-indexed)
  #     y.starts   : y-direction tile origin positions (standardized, 1-indexed)
  #   coords    : original grid coordinates
  #   pt.coords : (if pt.data provided) original point process coordinates

  ### 1. Grid setup (mirrors spbbGrid)

  orig.coords <- data[, coord.names]

  x.coords <- unique(orig.coords[, 1])
  delta.x  <- abs(base::diff(x.coords[1:2]))
  lims.x   <- c(min(x.coords), max(x.coords))
  n.x      <- length(x.coords)

  y.coords <- unique(orig.coords[, 2])
  delta.y  <- abs(base::diff(y.coords[1:2]))
  lims.y   <- c(min(y.coords), max(y.coords))
  n.y      <- length(y.coords)

  # Standardize: grid cell (j, i) has standardized coords (j, i), j in 1:n.x, i in 1:n.y
  new.coords <- orig.coords
  new.coords[, 1] <- (new.coords[, 1] - lims.x[1] + delta.x) / delta.x
  new.coords[, 2] <- (new.coords[, 2] - lims.y[1] + delta.y) / delta.y

  # grid.index[i, j] = original row index of the point at (x=j, y=i)
  coord.index <- order(new.coords[, 1], new.coords[, 2], decreasing = c(FALSE, FALSE))
  grid.index  <- matrix(coord.index, nrow = n.y, ncol = n.x)

  ### 2. Set up for point process data (optional)

  has.pt <- !is.null(pt.data) && !is.null(pt.coord.names)

  if (has.pt) {
    pt.orig.coords <- pt.data[, pt.coord.names, drop = FALSE]
    # Same standardization as the grid: integer k <-> grid cell k
    pt.std.x <- (pt.orig.coords[, 1] - lims.x[1] + delta.x) / delta.x
    pt.std.y <- (pt.orig.coords[, 2] - lims.y[1] + delta.y) / delta.y
  }

  ### 3. Build tiles for each block length

  tiles <- vector("list", length(block.l))

  for (m in seq_along(block.l)) {
    b.l <- block.l[m]

    if (b.l > n.x || b.l > n.y) {
      stop(sprintf(
        "block.l = %d exceeds grid dimensions (%d x %d).", b.l, n.x, n.y
      ))
    }

    # Step size between successive tile origins (in grid cells).
    # kappa = 0  -> step = b.l (non-overlapping)
    # kappa -> 1 -> step = 1  (maximally overlapping)
    step <- max(1L, round(b.l * (1 - kappa)))

    # Valid tile origin positions: tile must fit entirely within the grid.
    # If the grid dimension is not divisible by step, the last tile is shifted
    # back to end exactly at the boundary, ensuring full coverage.
    x.starts <- seq(1L, n.x - b.l + 1L, by = step)
    if (tail(x.starts, 1L) != n.x - b.l + 1L) x.starts <- c(x.starts, n.x - b.l + 1L)

    y.starts <- seq(1L, n.y - b.l + 1L, by = step)
    if (tail(y.starts, 1L) != n.y - b.l + 1L) y.starts <- c(y.starts, n.y - b.l + 1L)

    n.tiles.x  <- length(x.starts)
    n.tiles.y  <- length(y.starts)
    n.tiles    <- n.tiles.x * n.tiles.y
    n.per.tile <- b.l * b.l

    # Each column of tile.indices holds the original row indices for one tile.
    # Tiles are ordered: y-index varies fastest (column-major over tile origins).
    tile.indices <- matrix(NA_integer_, nrow = n.per.tile, ncol = n.tiles)

    if (has.pt) pt.indices <- vector("list", n.tiles)

    tile.num <- 1L
    for (j in seq_along(x.starts)) {
      for (i in seq_along(y.starts)) {
        x.block <- x.starts[j] + seq_len(b.l) - 1L   # column indices in grid.index
        y.block <- y.starts[i] + seq_len(b.l) - 1L   # row    indices in grid.index

        # grid.index[y.block, x.block] is b.l x b.l; c() flattens column-major
        tile.indices[, tile.num] <- c(grid.index[y.block, x.block])

        # Point process: collect indices of points within this tile's spatial extent.
        # Tile covers standardized x in [x.starts[j] - 0.5, x.starts[j] + b.l - 0.5)
        # and standardized y in [y.starts[i] - 0.5, y.starts[i] + b.l - 0.5).
        if (has.pt) {
          pt.indices[[tile.num]] <- which(
            pt.std.x >= x.starts[j] - 0.5 & pt.std.x < x.starts[j] + b.l - 0.5 &
            pt.std.y >= y.starts[i] - 0.5 & pt.std.y < y.starts[i] + b.l - 0.5
          )
        }

        tile.num <- tile.num + 1L
      }
    }

    out <- list(
      indices   = tile.indices,
      n.tiles   = n.tiles,
      n.tiles.x = n.tiles.x,
      n.tiles.y = n.tiles.y,
      step      = step,
      x.starts  = x.starts,
      y.starts  = y.starts
    )
    if (has.pt) out$pt.indices <- pt.indices
    tiles[[m]] <- out
  }

  names(tiles) <- paste0("bl.", block.l)

  result <- list(
    tiles  = tiles,
    coords = orig.coords
  )
  if (has.pt) result$pt.coords <- pt.orig.coords
  result
}
