### spbbGrid.R
### Block bootstrap for gridded raster data, with optional point process support.

spbbGrid <- function(
  data, coord.names = NULL, n.boot, block.l,
  pt.data = NULL, pt.coord.names = NULL, shift = TRUE
) {
  # Create a block bootstrap sample for a given lattice.
  # Optionally, simultaneously bootstraps a set of point process observations
  # (not on the grid) using the same block structure.
  #
  # Args:
  #   data          : gridded observations, as a terra SpatRaster, raster
  #                   Raster* object, sf object (POINT geometry, one point per
  #                   cell), sp Spatial* object, or data.frame/matrix. The grid
  #                   must be complete and regularly spaced (one observation per
  #                   cell); the number of x and y cells may differ.
  #   coord.names   : character(2), x and y column names; required only when
  #                   `data` is a data.frame/matrix
  #   n.boot        : number of bootstrap replicates
  #   block.l       : integer vector of block side lengths (in grid cells) to use
  #   shift         : if TRUE (default), the tiling is moved by a random whole
  #                   number of cells (0 to block.l - 1 in x and y) for each
  #                   replicate, so tile edges do not fall in the same places
  #                   every time; if FALSE, tiles start at the grid's corner
  #                   (moved at random only when block.l does not divide the
  #                   grid dimensions)
  #   pt.data       : (optional) point process observations, as an sf object
  #                   (POINT geometry), sp Spatial* object, or data.frame/matrix
  #   pt.coord.names: character(2), x and y column names; required only when
  #                   `pt.data` is a data.frame/matrix
  #
  # Returns a list with:
  #   bootstraps   : list (one per block length) of n.total x n.boot integer matrices
  #                  of indices into `data` (row indices for a data.frame/sf/sp
  #                  object, cell numbers for a raster)
  #   coords       : original grid coordinates
  #   pt.bootstraps: (if pt.data provided) list (one per block length) of length-n.boot
  #                  lists, each element a variable-length integer vector of row indices
  #                  into `pt.data`
  #   pt.coords    : (if pt.data provided) original point process coordinates

  ### 1. Set up for gridded data

  checkLonLat(list(data = data, pt.data = pt.data))
  orig.coords <- spCoords(data, coord.names, arg = "data")

  x.coords <- sort(unique(orig.coords[, 1]))
  y.coords <- sort(unique(orig.coords[, 2]))
  n.x      <- length(x.coords)
  n.y      <- length(y.coords)
  n.total  <- n.x * n.y

  if (n.x < 2 || n.y < 2) {
    stop("The grid must have at least two distinct x and y coordinates.", call. = FALSE)
  }
  if (nrow(orig.coords) != n.total || anyDuplicated(orig.coords) > 0) {
    stop(sprintf(
      "`data` must be a complete grid with one observation per cell: found %d observations for a %d x %d grid.",
      nrow(orig.coords), n.x, n.y
    ), call. = FALSE)
  }

  delta.x <- base::diff(x.coords[1:2])
  lims.x  <- range(x.coords)
  delta.y <- base::diff(y.coords[1:2])
  lims.y  <- range(y.coords)

  # Standardize grid coordinates so each unit equals one grid cell;
  # standardized x in {1, ..., n.x}, standardized y in {1, ..., n.y}
  new.coords <- orig.coords
  new.coords[, 1] <- (new.coords[, 1] - lims.x[1] + delta.x) / delta.x
  new.coords[, 2] <- (new.coords[, 2] - lims.y[1] + delta.y) / delta.y

  # Build an (n.y x n.x) matrix mapping (row = y_idx, col = x_idx) -> original row index
  coord.index <- order(new.coords[, 1], new.coords[, 2])
  grid.index  <- matrix(coord.index, nrow = n.y, ncol = n.x)

  ### 2. Set up for point process data (optional)

  has.pt <- !is.null(pt.data)

  if (has.pt) {
    pt.orig.coords <- spCoords(pt.data, pt.coord.names, arg = "pt.data")

    # Apply the same standardization as the grid so that
    # standardized coordinate k corresponds to grid cell k
    pt.std.x <- (pt.orig.coords[, 1] - lims.x[1] + delta.x) / delta.x
    pt.std.y <- (pt.orig.coords[, 2] - lims.y[1] + delta.y) / delta.y
    # Points outside [0.5, n.x+0.5] x [0.5, n.y+0.5] lie outside the grid domain
    # and will not appear in any bootstrap draw.
  }

  ### 3. Generate bootstrap samples

  if (missing(block.l) || length(block.l) == 0) {
    stop("`block.l` must be supplied.", call. = FALSE)
  }
  if (any(block.l < 1) || any(block.l > min(n.x, n.y))) {
    stop(sprintf(
      "`block.l` must be between 1 and %d (the smaller grid dimension).", min(n.x, n.y)
    ), call. = FALSE)
  }
  block.l <- as.integer(block.l)

  bootstraps <- vector("list", length(block.l))
  if (has.pt) pt.bootstraps <- vector("list", length(block.l))

  for (m in seq_along(block.l)) {

    b.l <- block.l[m]

    # one extra tile in each direction when the tiling is shifted
    n.tiles.x <- ceiling(n.x / b.l) + shift
    n.tiles.y <- ceiling(n.y / b.l) + shift
    n.blocks  <- n.tiles.x * n.tiles.y

    # Grid bootstrap: fixed n.total x n.boot matrix of indices
    boot.samples <- matrix(NA_integer_, nrow = n.total, ncol = n.boot)

    # Point process bootstrap: variable-length index vector per draw
    if (has.pt) pt.boot.samples <- vector("list", n.boot)

    for (k in 1:n.boot) {

      # Draw source-block lower-left corners for each tile position.
      # Column 1: row offset into grid.index (y direction)
      # Column 2: col offset into grid.index (x direction)
      sample.corners <- cbind(
        sample(n.y - b.l + 1, n.blocks, replace = TRUE),
        sample(n.x - b.l + 1, n.blocks, replace = TRUE)
      )

      # Build the over-sized "tiled" sample matrix (rows = y, cols = x) before trimming
      big.sample <- matrix(
        NA_integer_,
        nrow = n.tiles.y * b.l,
        ncol = n.tiles.x * b.l
      )

      for (i in 1:n.tiles.y) {
        for (j in 1:n.tiles.x) {

          y.block <- (i - 1) * b.l + 1:b.l   # row indices in big.sample (y dir)
          x.block <- (j - 1) * b.l + 1:b.l   # col indices in big.sample (x dir)

          bt.ij <- sample.corners[(i - 1) * n.tiles.x + j, ]

          big.sample[y.block, x.block] <- grid.index[
            bt.ij[1]:(bt.ij[1] + b.l - 1),
            bt.ij[2]:(bt.ij[2] + b.l - 1)
          ]
        }
      }

      # Randomly trim big.sample down to exactly n.y x n.x when b.l doesn't divide evenly
      # (with shift = TRUE, a uniform offset of 0 to b.l - 1 cells in each direction)
      start.pts <- c(1L, 1L)
      if (shift) {
        start.pts <- c(sample(b.l, 1), sample(b.l, 1))
      } else {
        if (n.y %% b.l > 0) start.pts[1] <- sample(nrow(big.sample) - n.y + 1, 1)
        if (n.x %% b.l > 0) start.pts[2] <- sample(ncol(big.sample) - n.x + 1, 1)
      }

      # Column-major flattening (y fastest, then x) matches the ordering of coord.index
      boot.samples[, k] <- c(big.sample[
        start.pts[1]:(start.pts[1] + n.y - 1),
        start.pts[2]:(start.pts[2] + n.x - 1)
      ])

      ### Point process bootstrap -------------------------------------------
      # For each tile that overlaps the trimmed window, identify the source-block
      # region in standardized coordinates and collect point process indices.
      if (has.pt) {

        pt.boot.k <- integer(0)

        for (i_tile in 1:n.tiles.y) {
          for (j_tile in 1:n.tiles.x) {

            # Tile's extent in big.sample rows (y dir) and cols (x dir)
            r_tile_start <- (i_tile - 1) * b.l + 1L
            r_tile_end   <-  i_tile      * b.l
            c_tile_start <- (j_tile - 1) * b.l + 1L
            c_tile_end   <-  j_tile      * b.l

            # Intersection of this tile with the trimmed window
            r_min <- max(r_tile_start, start.pts[1])
            r_max <- min(r_tile_end,   start.pts[1] + n.y - 1L)
            c_min <- max(c_tile_start, start.pts[2])
            c_max <- min(c_tile_end,   start.pts[2] + n.x - 1L)

            # Skip tiles entirely outside the used window
            if (r_min > r_max || c_min > c_max) next

            # Source block corner for this tile
            bt.ij <- sample.corners[(i_tile - 1L) * n.tiles.x + j_tile, ]

            # Map the overlap region to its position within the source block.
            # bt.ij[1] is a row offset in grid.index (y direction).
            # bt.ij[2] is a col offset in grid.index (x direction).
            src.row.min <- bt.ij[1] + (r_min - r_tile_start)
            src.row.max <- bt.ij[1] + (r_max - r_tile_start)
            src.col.min <- bt.ij[2] + (c_min - c_tile_start)
            src.col.max <- bt.ij[2] + (c_max - c_tile_start)

            # A point process point at standardized (std.x, std.y) belongs to grid
            # cell (col=round(std.x), row=round(std.y)).  We use half-cell margins
            # so cell k "owns" the interval [k - 0.5, k + 0.5).
            in.block <- which(
              pt.std.x >= src.col.min - 0.5 & pt.std.x < src.col.max + 0.5 &
              pt.std.y >= src.row.min - 0.5 & pt.std.y < src.row.max + 0.5
            )

            pt.boot.k <- c(pt.boot.k, in.block)
          }
        }

        pt.boot.samples[[k]] <- pt.boot.k
      }
      ### End point process bootstrap ----------------------------------------
    }

    bootstraps[[m]] <- boot.samples[order(coord.index), , drop = FALSE]
    if (has.pt) pt.bootstraps[[m]] <- pt.boot.samples
  }

  ### 4. Return results

  names(bootstraps) <- paste0("bl.", block.l)

  result <- list(
    bootstraps = bootstraps,
    coords     = orig.coords
  )

  if (has.pt) {
    names(pt.bootstraps) <- paste0("bl.", block.l)
    result$pt.bootstraps <- pt.bootstraps
    result$pt.coords     <- pt.orig.coords
  }

  result
}
