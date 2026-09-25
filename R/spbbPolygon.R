### spbbPolygon.R
### Block bootstrap for point data observed in an irregular (polygon) study area.

spbbPolygon <- function(
  data, boundary = NULL, n.boot, block.l, coord.names = NULL, shift = TRUE,
  output = c("objects", "indices")
){
  # Create block bootstrap samples of point data observed inside a polygon
  # (e.g. a state outline). The bounding box of the polygon is tiled with
  # square tiles of side `block.l`; every tile that overlaps the polygon is
  # filled with the points of a square source block drawn uniformly at random
  # from all squares lying entirely inside the polygon, translated onto the
  # tile. Points that land outside the polygon (in tiles crossing the
  # boundary) are dropped.
  #
  # Args:
  #   data       : point observations, as an sf/sfc object (POINT geometry),
  #                sp SpatialPoints* object, or data.frame/matrix
  #   boundary   : study area, as an sf/sfc (MULTI)POLYGON or sp Spatial
  #                polygons object; multiple features are unioned. If NULL
  #                (default), the convex hull of the points is used.
  #   n.boot     : number of bootstrap replicates
  #   block.l    : vector of block side lengths, in coordinate units
  #   coord.names: character(2), x and y column names; required only when
  #                `data` is a data.frame/matrix
  #   shift      : if TRUE (default), the tiling is moved by a random offset
  #                (uniform in [0, block.l) in x and y) for each replicate, so
  #                tile edges do not fall in the same places every time; if
  #                FALSE, tiles always start at the bounding box's lower-left
  #                corner
  #   output     : "objects" (default) returns each replicate as an object of
  #                the same class as `data`; "indices" returns each replicate
  #                in compact form, list(idx = rows of `data`, xy = new
  #                coordinates), which is far smaller (e.g. for saving to disk
  #                and running replicates as separate cluster jobs). Turn a
  #                compact replicate into a full object with buildReplicate().
  #
  # Returns a named list (one element per block length), each a length-n.boot
  # list of bootstrap samples with the same class as `data`, containing copies
  # of the resampled observations moved to their new locations.
  # With output = "indices", each sample is instead list(idx, xy).

  output <- match.arg(output)
  checkPointData(data)
  checkLonLat(list(data = data, boundary = boundary))
  reps <- polygonReplicates(data, boundary, n.boot, block.l, coord.names, shift)
  if (output == "indices") return(reps)
  lapply(reps, lapply, function(r) movePoints(data, r$idx, r$xy, coord.names))
}


polygonReplicates <- function(
  data, boundary = NULL, n.boot, block.l, coord.names = NULL, shift = TRUE
){
  # Workhorse for spbbPolygon(): returns, for each block length, a length-n.boot
  # list of replicates in compact form, list(idx = rows of `data`, xy = new
  # coordinates), so that callers can build full objects one at a time.

  checkPointData(data)
  pts <- as.matrix(spCoords(data, coord.names, arg = "data"))

  boundary <- resolveBoundary(boundary, data, pts)
  bb <- sf::st_bbox(boundary)

  if (missing(block.l) || length(block.l) == 0) {
    stop("`block.l` must be supplied.", call. = FALSE)
  }
  max.l <- min(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  if (!is.numeric(block.l) || anyNA(block.l) || any(block.l <= 0) || any(block.l > max.l)) {
    stop(sprintf(
      "`block.l` must be greater than 0 and at most %g (the shorter side of the boundary's bounding box).",
      max.l
    ), call. = FALSE)
  }

  n.out <- sum(!inPolygon(pts, boundary))
  if (n.out > 0) {
    warning(sprintf(
      "%d of %d points lie outside `boundary` and will never be resampled.",
      n.out, nrow(pts)
    ), call. = FALSE)
  }

  # points sorted by x, so the points in any square are found by binary search
  x.ord <- order(pts[, 1])
  xs    <- pts[x.ord, 1]
  ys    <- pts[x.ord, 2]

  bootstraps <- vector("list", length(block.l))

  for (m in seq_along(block.l)) {
    b.l <- block.l[m]

    ### 1. Destination tiles for each replicate, keeping those overlapping the
    ###    boundary and noting which lie entirely inside it (their points never
    ###    need clipping). With shift = TRUE the tiling is moved by a random
    ###    offset in [0, b.l) in x and y for each replicate.
    n.tiles.x <- ceiling((bb[["xmax"]] - bb[["xmin"]]) / b.l) + shift
    n.tiles.y <- ceiling((bb[["ymax"]] - bb[["ymin"]]) / b.l) + shift
    checkTileCount(n.tiles.x * n.tiles.y, b.l)
    grid.x <- (rep(seq_len(n.tiles.x), times = n.tiles.y) - 1) * b.l
    grid.y <- (rep(seq_len(n.tiles.y), each = n.tiles.x) - 1) * b.l

    off.x <- off.y <- numeric(n.boot)
    if (shift) {
      off.x <- stats::runif(n.boot, 0, b.l)
      off.y <- stats::runif(n.boot, 0, b.l)
    }
    layouts <- vector("list", n.boot)
    for (k in seq_len(n.boot)) {
      layouts[[k]] <- if (k > 1 && !shift) layouts[[1]] else tileLayout(
        bb[["xmin"]] - off.x[k] + grid.x, bb[["ymin"]] - off.y[k] + grid.y, b.l, boundary
      )
    }
    n.active <- vapply(layouts, function(l) length(l$x), 1L)

    ### 2. Source blocks: one square inside the boundary per tile per replicate
    corners <- sampleSquares(boundary, b.l, sum(n.active))

    # index range of points with src.x <= x < src.x + b.l (y is checked below)
    lo <- findInterval(corners[, 1], xs, left.open = TRUE) + 1L
    hi <- findInterval(corners[, 1] + b.l, xs, left.open = TRUE)

    ### 3. Assemble each replicate, dropping points moved outside the boundary
    ###    (only points in tiles crossing the boundary need testing)
    reps <- vector("list", n.boot)
    d <- 0L
    for (k in seq_len(n.boot)) {
      lay  <- layouts[[k]]
      rows <- d + seq_len(n.active[k])     # this replicate's rows of `corners`
      idx.k <- vector("list", n.active[k])
      for (t in seq_len(n.active[k])) {
        d <- d + 1L
        if (hi[d] < lo[d]) next
        cand <- lo[d]:hi[d]
        idx.k[[t]] <- cand[ys[cand] >= corners[d, 2] & ys[cand] < corners[d, 2] + b.l]
      }
      n.k <- lengths(idx.k)
      j   <- unlist(idx.k, use.names = FALSE)
      tt  <- rep(seq_len(n.active[k]), n.k)
      xy  <- cbind(
        xs[j] - corners[rows[tt], 1] + lay$x[tt],
        ys[j] - corners[rows[tt], 2] + lay$y[tt]
      )
      keep <- rep(TRUE, length(j))
      on.edge <- lay$edge[tt]
      if (any(on.edge)) keep[on.edge] <- inPolygon(xy[on.edge, , drop = FALSE], boundary)
      reps[[k]] <- list(idx = x.ord[j[keep]], xy = xy[keep, , drop = FALSE])
    }
    bootstraps[[m]] <- reps
  }

  names(bootstraps) <- paste0("bl.", block.l)
  bootstraps
}


resolveBoundary <- function(boundary, data, pts) {
  # Return the study area as a single-feature polygon sfc in the coordinates of
  # `data`. The block bootstrap works in planar coordinates (squares and
  # shifts in coordinate units, even for longitude/latitude), so the CRS is
  # dropped after any transformation: all geometry is then done by GEOS in
  # the plane, never by s2 on the sphere, and areas are in coordinate units.
  crs <- if (inherits(data, c("sf", "sfc"))) {
    sf::st_crs(data)
  } else if (inherits(data, "Spatial")) {
    sf::st_crs(methods::slot(data, "proj4string"))
  } else {
    sf::NA_crs_
  }

  if (is.null(boundary)) {
    warning(
      "No `boundary` supplied; using the convex hull of the points. Supplying the study area is best when it is known.",
      call. = FALSE
    )
    hull <- pts[grDevices::chull(pts), , drop = FALSE]
    return(sf::st_sfc(sf::st_polygon(list(rbind(hull, hull[1, ])))))
  }

  if (inherits(boundary, "Spatial")) boundary <- sf::st_as_sf(boundary)
  if (!inherits(boundary, c("sf", "sfc"))) {
    stop("`boundary` must be an sf/sfc polygon or an sp Spatial polygons object.", call. = FALSE)
  }
  geom <- sf::st_geometry(boundary)
  if (!all(as.character(sf::st_geometry_type(geom)) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`boundary` must have POLYGON or MULTIPOLYGON geometry.", call. = FALSE)
  }

  if (!is.na(crs) && !is.na(sf::st_crs(geom)) && sf::st_crs(geom) != crs) {
    geom <- sf::st_transform(geom, crs)
  }
  geom <- sf::st_set_crs(geom, NA)
  if (length(geom) > 1) geom <- sf::st_union(geom)
  geom
}


tileLayout <- function(tile.x, tile.y, b.l, boundary) {
  # Of the tiles with lower-left corners (tile.x, tile.y), keep those that
  # overlap `boundary`; `edge` flags tiles not lying entirely inside it.
  tiles    <- makeSquares(tile.x, tile.y, b.l, sf::st_crs(boundary))
  active   <- seq_along(tiles) %in% sf::st_intersects(boundary, tiles)[[1]]
  interior <- seq_along(tiles) %in% sf::st_contains(boundary, tiles)[[1]]
  list(x = tile.x[active], y = tile.y[active], edge = !interior[active])
}


makeSquares <- function(x, y, b.l, crs) {
  # Axis-aligned squares with lower-left corners (x, y) and side b.l, as sfc.
  sf::st_sfc(lapply(seq_along(x), function(i) {
    sf::st_polygon(list(matrix(
      c(x[i], x[i] + b.l, x[i] + b.l, x[i], x[i],
        y[i], y[i], y[i] + b.l, y[i] + b.l, y[i]),
      ncol = 2
    )))
  }), crs = crs)
}


inPolygon <- function(xy, boundary) {
  # Logical: which rows of the coordinate matrix `xy` fall inside `boundary`.
  if (nrow(xy) == 0) return(logical(0))
  p <- sf::st_geometry(sf::st_as_sf(as.data.frame(xy), coords = 1:2, crs = sf::st_crs(boundary)))
  # boundary first, so GEOS prepares the polygon once for all points
  seq_len(nrow(xy)) %in% sf::st_intersects(boundary, p)[[1]]
}


validCorners <- function(boundary, b.l) {
  # Region of lower-left corners c such that the square [c, c + b.l]^2 lies
  # entirely inside `boundary`. The square is inside exactly when its corner
  # is inside and the square does not touch the boundary; the corners whose
  # square touches boundary segment [p, q] form the convex hull of
  # {p, q} - [0, b.l]^2 (a hexagon), so the valid region is the boundary minus
  # the union of these hexagons, computed once per block length.
  xy  <- sf::st_coordinates(sf::st_boundary(boundary))
  grp <- do.call(paste, as.data.frame(xy[, -(1:2), drop = FALSE]))
  nxt <- c(seq_len(nrow(xy))[-1], NA)
  seg <- which(!is.na(nxt) & grp == c(grp[-1], NA))
  # skip zero-length segments (repeated vertices), which add nothing
  seg <- seg[xy[seg, 1] != xy[nxt[seg], 1] | xy[seg, 2] != xy[nxt[seg], 2]]

  off <- cbind(c(0, -b.l, -b.l, 0), c(0, 0, -b.l, -b.l))
  hexagons <- lapply(seg, function(i) {
    v <- rbind(sweep(off, 2, xy[i, 1:2], `+`), sweep(off, 2, xy[nxt[i], 1:2], `+`))
    h <- v[grDevices::chull(v), , drop = FALSE]
    sf::st_polygon(list(rbind(h, h[1, ])))
  })
  touching <- sf::st_union(sf::st_sfc(hexagons, crs = sf::st_crs(boundary)))
  sf::st_difference(boundary, touching)
}


sampleSquares <- function(boundary, b.l, n) {
  # Lower-left corners of `n` squares of side b.l drawn uniformly from all
  # squares lying entirely inside `boundary`: points drawn uniformly over the
  # region of valid corners (validCorners) by rejection sampling.
  valid <- validCorners(boundary, b.l)
  if (length(valid) == 0 || sf::st_is_empty(valid) || as.numeric(sf::st_area(valid)) <= 0) {
    stop(sprintf(
      "No square of side %g fits inside `boundary`; use a smaller block length.", b.l
    ), call. = FALSE)
  }
  vb <- sf::st_bbox(valid)
  rate <- as.numeric(sf::st_area(valid)) /
    ((vb[["xmax"]] - vb[["xmin"]]) * (vb[["ymax"]] - vb[["ymin"]]))

  out <- matrix(numeric(0), ncol = 2)
  while (nrow(out) < n) {
    m  <- ceiling(1.1 * (n - nrow(out)) / rate) + 10
    xy <- cbind(stats::runif(m, vb[["xmin"]], vb[["xmax"]]),
                stats::runif(m, vb[["ymin"]], vb[["ymax"]]))
    out <- rbind(out, xy[inPolygon(xy, valid), , drop = FALSE])
  }
  out[seq_len(n), , drop = FALSE]
}
