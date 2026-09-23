### spbbGrid.R
### Block bootstrap for gridded raster data, with optional point process support.

spbbGrid <- function(
  data, coord.names, n.boot, block.l = NULL,
  pt.data = NULL, pt.coord.names = NULL
) {
  # Create a block bootstrap sample for a given lattice.
  # Optionally, simultaneously bootstraps a set of point process observations
  # (not on the grid) using the same block structure.
  #
  # Args:
  #   data          : data.frame with gridded observations
  #   coord.names   : character(2), column names for x and y coordinates in `data`
  #   n.boot        : number of bootstrap replicates
  #   block.l       : integer vector of block side lengths to use
  #   pt.data       : (optional) data.frame of point process observations
  #   pt.coord.names: (optional) character(2), column names for x and y in `pt.data`
  #
  # Returns a list with:
  #   bootstraps   : list (one per block length) of n.total x n.boot integer matrices
  #                  of row indices into `data`
  #   coords       : original grid coordinates
  #   pt.bootstraps: (if pt.data provided) list (one per block length) of length-n.boot
  #                  lists, each element a variable-length integer vector of row indices
  #                  into `pt.data`
  #   pt.coords    : (if pt.data provided) original point process coordinates

  ### 1. Set up for gridded data

  orig.coords <- data[, coord.names]

  x.coords <- unique(orig.coords[, 1])
  delta.x  <- abs(base::diff(x.coords[1:2]))
  lims.x   <- c(min(x.coords), max(x.coords))
  n.x      <- length(x.coords)

  y.coords <- unique(orig.coords[, 2])
  delta.y  <- abs(base::diff(y.coords[1:2]))
  lims.y   <- c(min(y.coords), max(y.coords))
  n.y      <- length(y.coords)

  # Standardize grid coordinates so each unit equals one grid cell;
  # standardized x in {1, ..., n.x}, standardized y in {1, ..., n.y}
  new.coords <- orig.coords
  new.coords[, 1] <- (new.coords[, 1] - lims.x[1] + delta.x) / delta.x
  new.coords[, 2] <- (new.coords[, 2] - lims.y[1] + delta.y) / delta.y

  n.total <- n.x * n.y

  # Build an (n.y x n.x) matrix mapping (row=y_idx, col=x_idx) -> original row index
  coord.index <- order(new.coords[, 1], new.coords[, 2], decreasing = c(FALSE, FALSE))
  grid.index  <- matrix(coord.index, nrow = n.y, ncol = n.x)

  ### 2. Set up for point process data (optional)

  has.pt <- !is.null(pt.data) && !is.null(pt.coord.names)

  if (has.pt) {
    pt.orig.coords <- pt.data[, pt.coord.names, drop = FALSE]

    # Apply the same standardization as the grid so that
    # standardized coordinate k corresponds to grid cell k
    pt.std.x <- (pt.orig.coords[, 1] - lims.x[1] + delta.x) / delta.x
    pt.std.y <- (pt.orig.coords[, 2] - lims.y[1] + delta.y) / delta.y
    # Points outside [0.5, n.x+0.5] x [0.5, n.y+0.5] lie outside the grid domain
    # and will not appear in any bootstrap draw.
  }

  ### 3. Generate bootstrap samples

  bootstraps <- vector("list", length(block.l))
  if (has.pt) pt.bootstraps <- vector("list", length(block.l))

  for (m in seq_along(block.l)) {

    b.l <- block.l[m]

    n.tiles.x <- ceiling(n.x / b.l)
    n.tiles.y <- ceiling(n.y / b.l)
    n.blocks  <- n.tiles.x * n.tiles.y

    # Grid bootstrap: fixed n.total x n.boot matrix of indices
    boot.samples <- matrix(NA_integer_, nrow = n.total, ncol = n.boot)

    # Point process bootstrap: variable-length index vector per draw
    if (has.pt) pt.boot.samples <- vector("list", n.boot)

    for (k in 1:n.boot) {

      # Draw source-block lower-left corners for each tile position.
      # Column 1: row offset into grid.index  (bounded by n.x due to existing convention)
      # Column 2: col offset into grid.index  (bounded by n.y due to existing convention)
      sample.corners <- cbind(
        sample(n.x - b.l + 1, n.blocks, replace = TRUE),
        sample(n.y - b.l + 1, n.blocks, replace = TRUE)
      )

      # Build the over-sized "tiled" sample matrix before trimming
      big.sample <- matrix(
        NA_integer_,
        nrow = n.tiles.x * b.l,
        ncol = n.tiles.y * b.l
      )

      for (i in 1:n.tiles.y) {
        for (j in 1:n.tiles.x) {

          y.block <- (i - 1) * b.l + 1:b.l   # column indices in big.sample (y dir)
          x.block <- (j - 1) * b.l + 1:b.l   # row    indices in big.sample (x dir)

          bt.ij <- sample.corners[(i - 1) * n.tiles.x + j, ]

          big.sample[x.block, y.block] <- grid.index[
            bt.ij[1]:(bt.ij[1] + b.l - 1),
            bt.ij[2]:(bt.ij[2] + b.l - 1)
          ]
        }
      }

      # Randomly trim big.sample down to exactly n.x x n.y when b.l doesn't divide evenly
      start.pts <- c(1L, 1L)
      if (n.x %% b.l > 0) start.pts[1] <- sample(nrow(big.sample) - n.x + 1, 1)
      if (n.y %% b.l > 0) start.pts[2] <- sample(ncol(big.sample) - n.y + 1, 1)

      boot.samples[, k] <- c(big.sample[
        start.pts[1]:(start.pts[1] + n.x - 1),
        start.pts[2]:(start.pts[2] + n.y - 1)
      ])

      ### Point process bootstrap -------------------------------------------
      # For each tile that overlaps the trimmed window, identify the source-block
      # region in standardized coordinates and collect point process indices.
      if (has.pt) {

        pt.boot.k <- integer(0)

        for (i_tile in 1:n.tiles.y) {
          for (j_tile in 1:n.tiles.x) {

            # Tile's extent in big.sample rows (x dir) and cols (y dir)
            r_tile_start <- (j_tile - 1) * b.l + 1L
            r_tile_end   <-  j_tile      * b.l
            c_tile_start <- (i_tile - 1) * b.l + 1L
            c_tile_end   <-  i_tile      * b.l

            # Intersection of this tile with the trimmed window
            r_min <- max(r_tile_start, start.pts[1])
            r_max <- min(r_tile_end,   start.pts[1] + n.x - 1L)
            c_min <- max(c_tile_start, start.pts[2])
            c_max <- min(c_tile_end,   start.pts[2] + n.y - 1L)

            # Skip tiles entirely outside the used window
            if (r_min > r_max || c_min > c_max) next

            # Source block corner for this tile
            bt.ij <- sample.corners[(i_tile - 1L) * n.tiles.x + j_tile, ]

            # Map the overlap region to its position within the source block.
            # bt.ij[1] is a row offset in grid.index (y direction in original code).
            # bt.ij[2] is a col offset in grid.index (x direction in original code).
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

    bootstraps[[m]] <- boot.samples[order(coord.index), ]
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


inBox <- function(pts, corner, b.l){
  which(pts[,1] >= corner[1] & pts[,1] <= corner[1] + b.l &
          pts[,2] >= corner[2] & pts[,2] <= corner[2] + b.l)
}


spbbGrid.pt <- function(
  data, coord.names, n.boot, block.l = NULL, box.x = c(0,100), box.y = c(0,100)
){

  # grab data, determine box size
  pts <- data[, coord.names]
  min.x <- min(box.x); max.x <- max(box.x); l.x <- max.x - min.x
  min.y <- min(box.y); max.y <- max(box.y); l.y <- max.y - min.y

  # save bootstrap results
  bootstraps <- vector("list", length(block.l))

  for (m in seq_along(block.l)) {
    # block length
    b.l <- as.integer(block.l[m])

    # number of blocks
    n.tiles.x <- ceiling(l.x / b.l)
    n.tiles.y <- ceiling(l.y / b.l)
    n.blocks  <- n.tiles.x * n.tiles.y

    # save bootstrap samples
    boot.samples <- list()

    for (k in 1:n.boot) {

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
          orig.pts <- pts[new.pts,]

          # standardize the points
          std.pts <- cbind(
            orig.pts$x - corners.k[curr.idx,1], orig.pts$y - corners.k[curr.idx,2]
          )

          # move points into new tiled grid
          pts.lst[[curr.idx]] <- cbind((i-1) * b.l + std.pts[,1], (j - 1) * b.l + std.pts[,2])
          pts.idx.lst[[curr.idx]] <- new.pts
        }
      }

      # convert list to data frame
      pts.full <- do.call(rbind, pts.lst)
      pts.idx <- unlist(pts.idx.lst)

      # remove any points outside the boundary
      keep.pts <- which(pts.full[,1] <= max(box.x) & pts.full[,2] <= max(box.y) &
                          pts.full[,1] >= min(box.x) & pts.full[,2] >= min(box.y))
      pts.k <- pts.full[keep.pts,]
      pts.idx.k <- pts.idx[keep.pts]

      df.k <- data[pts.idx.k, ]
      df.k[,coord.names] <- pts.k

      # save results
      boot.samples[[k]] <- df.k
    }
    bootstraps[[m]] <- boot.samples
  }

  names(bootstraps) <- paste0("bl.", block.l)
  return(bootstraps)
}

