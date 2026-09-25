### blockLength.R
### Empirical block length selection (Nordman & Lahiri style) for the spatial
### block bootstrap of a user-supplied statistic.

blockLength <- function(
  data, statistic, n.boot, type = c("grid", "point", "polygon"),
  coord.names = NULL, box.x = NULL, box.y = NULL, boundary = NULL,
  shift = TRUE, c1 = 1, c2 = 0.5, overlap = TRUE, apply.fun = lapply, ...
) {
  # Choose the block length for spboot() by the nonparametric plug-in method:
  # bootstrap variances at three pilot block lengths, b1, b2 and 2 * b2, give
  # estimates of the variance and bias of the block bootstrap variance
  # estimator, which are plugged into the MSE-optimal block length.
  #
  # The procedure is carried out in "count" units (grid cells, or for point
  # data the average spacing between points, sqrt(area / n), where the area is
  # that of the box or boundary polygon). Block lengths for point data are
  # converted to coordinate units at the end.
  #
  # Args:
  #   data, statistic, type, coord.names, box.x, box.y, boundary, shift,
  #   apply.fun, ...:
  #             as in spboot()
  #   n.boot  : number of bootstrap replicates at each pilot block length
  #   c1, c2  : pilot block constants, b1 = c1 * n^(1/4) and b2 = c2 * n^(1/6)
  #   overlap : if TRUE, apply the overlapping-block correction sqrt(3/2)
  #
  # Returns an object of class "spbootBlockLength": a list with
  #   block.l: selected block length for each element of the statistic, ready
  #            to pass to spboot() (integer grid cells for type = "grid",
  #            coordinate units for type = "point")
  #   pilot  : the three pilot block lengths (same units as block.l)
  #   var    : 3 x p matrix of bootstrap variances at the pilot block lengths
  #   B.hat  : estimated bias constant for each element
  #   n      : number of grid cells or points
  #   unit   : length of one count unit in coordinate units (1 for grids)
  #   fits   : the three pilot spboot objects

  type <- match.arg(type)
  checkLonLat(list(data = data, boundary = boundary))
  pts  <- spCoords(data, coord.names, arg = "data")
  n    <- nrow(pts)

  if (type == "grid") {
    unit  <- 1
    max.l <- min(length(unique(pts[, 1])), length(unique(pts[, 2])))
  } else if (type == "point") {
    box   <- resolveBox(pts, box.x, box.y)
    box.x <- box$x; box.y <- box$y
    area  <- diff(range(box.x)) * diff(range(box.y))
    unit  <- sqrt(area / n)
    max.l <- min(diff(range(box.x)), diff(range(box.y)))
  } else {
    boundary <- resolveBoundary(boundary, data, as.matrix(pts))
    bb    <- sf::st_bbox(boundary)
    area  <- as.numeric(sf::st_area(boundary))
    unit  <- sqrt(area / n)
    max.l <- min(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  }

  ### 1. Pilot block lengths (count units)

  b1 <- c1 * n^(1/4)
  b2 <- c2 * n^(1/6)
  if (type == "grid") {
    b1 <- max(1, round(b1))
    b2 <- max(1, round(b2))
  }
  pilot <- c(b1 = b1, b2 = b2, b3 = 2 * b2) * unit

  if (any(pilot > max.l)) {
    stop(sprintf(
      "Pilot block lengths (%s) exceed the largest allowed block length (%g); reduce `c1` or `c2`.",
      paste(signif(pilot, 4), collapse = ", "), max.l
    ), call. = FALSE)
  }

  ### 2. Bootstrap variances at each pilot block length

  fits <- lapply(pilot, function(b.l) {
    spboot(data, statistic, n.boot = n.boot, block.l = b.l, type = type,
           coord.names = coord.names, box.x = box.x, box.y = box.y,
           boundary = boundary, shift = shift, apply.fun = apply.fun, ...)
  })
  v <- do.call(rbind, lapply(fits, function(f) f$se^2))
  rownames(v) <- names(pilot)

  ### 3. Plug-in optimal block length

  # bias constant from the difference between block lengths b2 and 2 * b2
  B.hat <- 2 * b2 * (v["b2", ] - v["b3", ])

  # MSE-optimal length for non-overlapping blocks, then overlap correction
  b.opt <- (B.hat^2 * n / 2 / v["b1", ]^2)^(1/4)
  if (overlap) b.opt <- b.opt * sqrt(3 / 2)

  b.opt <- b.opt * unit
  if (type == "grid") b.opt <- ceiling(b.opt)

  # row extraction drops names when the statistic has a single element
  names(b.opt) <- names(B.hat) <- colnames(v)

  # keep the result usable: at least one count unit, at most the grid / box
  min.l   <- unit
  bad     <- !is.finite(b.opt)
  clamped <- !bad & (b.opt < min.l | b.opt > max.l)
  if (any(bad)) {
    warning(sprintf(
      "Block length could not be estimated for %s (non-finite variance estimates); using b2 = %g.",
      paste(names(b.opt)[bad], collapse = ", "), pilot[["b2"]]
    ), call. = FALSE)
    b.opt[bad] <- pilot[["b2"]]
  }
  if (any(clamped)) {
    warning(sprintf(
      "Estimated block length for %s fell outside [%g, %g] and was truncated.",
      paste(names(b.opt)[clamped], collapse = ", "), min.l, max.l
    ), call. = FALSE)
    b.opt <- pmin(pmax(b.opt, min.l), max.l)
  }

  structure(
    list(
      block.l = b.opt,
      pilot   = pilot,
      var     = v,
      B.hat   = B.hat,
      n       = n,
      unit    = unit,
      type    = type,
      fits    = fits
    ),
    class = "spbootBlockLength"
  )
}


print.spbootBlockLength <- function(x, ...) {
  cat(sprintf(
    "Block length selection (%s data, n = %d)\nPilot block lengths: %s\n\n",
    x$type, x$n, paste(signif(x$pilot, 4), collapse = ", ")
  ))
  print(data.frame(
    block.l = x$block.l,
    se.b1 = sqrt(x$var["b1", ]), se.b2 = sqrt(x$var["b2", ]), se.b3 = sqrt(x$var["b3", ]),
    row.names = names(x$block.l)
  ), ...)
  invisible(x)
}
