### spboot.R
### Spatial block bootstrap of a user-supplied statistic.

spboot <- function(
  data, statistic, n.boot, block.l, type = c("grid", "point", "polygon"),
  coord.names = NULL, box.x = NULL, box.y = NULL, boundary = NULL,
  shift = TRUE, apply.fun = lapply, ...
) {
  # Apply `statistic` to spatial block bootstrap replicates of `data` and
  # estimate its standard error.
  #
  # Each replicate is passed to `statistic` in the same form as `data`, with
  # every resampled observation relocated to its new position: a block drawn
  # from the upper right and placed in the lower left carries lower-left
  # coordinates. Spatial terms in the statistic (e.g. s(x, y)) therefore see
  # the bootstrap geometry, not the original one.
  #
  # Args:
  #   data       : for type = "grid", anything spbbGrid() accepts (SpatRaster,
  #                Raster*, sf, sp, data.frame/matrix); for type = "point" or
  #                "polygon", point data (sf, sp, data.frame/matrix)
  #   statistic  : function(data, ...) returning a numeric vector (e.g. an
  #                ATE, or regression coefficients); names are kept
  #   n.boot     : number of bootstrap replicates
  #   block.l    : block side length, either a single value used for every
  #                element of the statistic, or one value per element (named
  #                to match the statistic's names, or in the same order).
  #                Grid cells for type = "grid", coordinate units otherwise.
  #   type       : "grid" (spbbGrid), "point" (spbbPoint: rectangular study
  #                area) or "polygon" (spbbPolygon: irregular study area)
  #   coord.names: character(2), x and y column names; required only when
  #                `data` is a data.frame/matrix
  #   box.x/box.y: observation window for type = "point" (see spbbPoint())
  #   boundary   : study area polygon for type = "polygon" (see spbbPolygon())
  #   shift      : if TRUE (default), the tiling is randomly offset for each
  #                replicate; if FALSE, tile edges are the same in every
  #                replicate (see spbbGrid(), spbbPoint(), spbbPolygon())
  #   apply.fun  : lapply-like function used to evaluate the replicates, e.g.
  #                function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 4).
  #                All resampling happens before apply.fun is called, so the
  #                replicates do not depend on how they are evaluated. Point
  #                replicates are stored compactly and each is built only when
  #                `statistic` is applied to it.
  #   ...        : further arguments passed to `statistic`
  #
  # Returns an object of class "spboot": a list with
  #   t0     : statistic evaluated on `data`
  #   t      : n.boot x length(t0) matrix of replicate values
  #   se     : bootstrap standard errors (sd of each column of t)
  #   block.l: block length used for each element of t0
  #   n.boot, type
  #   n.error: number of replicates for which `statistic` threw an error
  #            (recorded as NA in t and ignored in se)

  type      <- match.arg(type)
  statistic <- match.fun(statistic)
  checkLonLat(list(data = data, boundary = boundary))

  # original estimate
  t0 <- checkStatistic(statistic(data, ...))
  p  <- length(t0)

  # one block length per element of the statistic
  block.l <- matchBlockLengths(block.l, names(t0))

  # settle the study area once (so any default-area warning is given once)
  if (type == "point") {
    box <- resolveBox(spCoords(data, coord.names, arg = "data"), box.x, box.y)
  } else if (type == "polygon") {
    boundary <- resolveBoundary(boundary, data,
                                as.matrix(spCoords(data, coord.names, arg = "data")))
  }

  t <- matrix(NA_real_, nrow = n.boot, ncol = p, dimnames = list(NULL, names(t0)))
  n.error <- stats::setNames(integer(p), names(t0))
  err.msg <- character(0)

  # run once per distinct block length, keeping the matching elements
  for (b.l in unique(block.l)) {
    cols <- which(block.l == b.l)

    if (type == "grid") {
      idx <- spbbGrid(data, coord.names, n.boot = n.boot, block.l = b.l,
                      shift = shift)$bootstraps[[1]]
      getReplicate <- function(k) gridReplicate(data, idx[, k], coord.names)
    } else {
      reps <- if (type == "point") {
        pointReplicates(data, coord.names, n.boot = n.boot, block.l = b.l,
                        box.x = box$x, box.y = box$y, shift = shift)[[1]]
      } else {
        polygonReplicates(data, boundary, n.boot = n.boot, block.l = b.l,
                          coord.names = coord.names, shift = shift)[[1]]
      }
      getReplicate <- function(k) movePoints(data, reps[[k]]$idx, reps[[k]]$xy, coord.names)
    }

    res <- apply.fun(seq_len(n.boot), function(k) {
      tryCatch(
        statistic(getReplicate(k), ...),
        error = function(e) structure(conditionMessage(e), class = "spbootError")
      )
    })

    for (k in seq_len(n.boot)) {
      if (inherits(res[[k]], "spbootError")) {
        n.error[cols] <- n.error[cols] + 1L
        err.msg <- c(err.msg, unclass(res[[k]]))
        next
      }
      t.k <- as.numeric(res[[k]])
      if (length(t.k) != p) {
        stop(sprintf(
          "`statistic` returned %d value(s) for a replicate but %d for the original data.",
          length(t.k), p
        ), call. = FALSE)
      }
      t[k, cols] <- t.k[cols]
    }
  }

  if (length(err.msg) > 0) {
    warning(sprintf(
      "`statistic` failed on %d replicate(s); these are NA in `t` and ignored in `se`. First error: %s",
      length(err.msg), err.msg[1]
    ), call. = FALSE)
  }

  structure(
    list(
      t0      = t0,
      t       = t,
      se      = apply(t, 2, stats::sd, na.rm = TRUE),
      block.l = block.l,
      n.boot  = n.boot,
      type    = type,
      shift   = shift,
      n.error = n.error
    ),
    class = "spboot"
  )
}


print.spboot <- function(x, ...) {
  cat(sprintf("Spatial block bootstrap (%s data, %d replicates)\n\n", x$type, x$n.boot))
  tab <- data.frame(
    estimate = x$t0, se = x$se, block.l = x$block.l,
    row.names = names(x$t0)
  )
  if (any(x$n.error > 0)) tab$failed <- x$n.error
  print(tab, ...)
  invisible(x)
}


checkStatistic <- function(t0) {
  # Validate the statistic's value on the original data and name its elements.
  if (!is.numeric(t0) || length(t0) == 0) {
    stop("`statistic` must return a non-empty numeric vector.", call. = FALSE)
  }
  t0 <- stats::setNames(as.numeric(t0), names(t0))
  if (is.null(names(t0))) {
    names(t0) <- if (length(t0) == 1) "t1" else paste0("t", seq_along(t0))
  }
  t0
}


matchBlockLengths <- function(block.l, stat.names) {
  # Expand `block.l` to one value per element of the statistic.
  p <- length(stat.names)
  if (missing(block.l) || length(block.l) == 0) {
    stop("`block.l` must be supplied.", call. = FALSE)
  }
  if (length(block.l) == 1) {
    block.l <- rep(unname(block.l), p)
  } else if (length(block.l) != p) {
    stop(sprintf(
      "`block.l` must have length 1 or %d (one per element of the statistic).", p
    ), call. = FALSE)
  } else if (!is.null(names(block.l))) {
    if (!setequal(names(block.l), stat.names)) {
      stop("Names of `block.l` must match the names of the statistic.", call. = FALSE)
    }
    block.l <- unname(block.l[stat.names])
  }
  stats::setNames(as.numeric(block.l), stat.names)
}


gridReplicate <- function(data, idx, coord.names = NULL) {
  # Build a grid bootstrap replicate: keep the geometry of `data` and move the
  # attribute values of rows (or cells) `idx` into positions 1, 2, ..., n.

  if (inherits(data, "SpatRaster")) {
    terra::setValues(data, terra::values(data)[idx, , drop = FALSE])

  } else if (inherits(data, "Raster")) {
    v <- raster::getValues(data)
    raster::setValues(data, if (is.matrix(v)) v[idx, , drop = FALSE] else v[idx])

  } else if (inherits(data, "sfc")) {
    stop("An sfc object has no attributes to resample; use an sf object.", call. = FALSE)

  } else if (inherits(data, "sf")) {
    out <- data[idx, ]
    sf::st_geometry(out) <- sf::st_geometry(data)
    rownames(out) <- NULL
    out

  } else if (inherits(data, "Spatial")) {
    if (!methods::.hasSlot(data, "data")) {
      stop("sp input must carry attributes (a Spatial*DataFrame).", call. = FALSE)
    }
    out <- data
    out@data <- data@data[idx, , drop = FALSE]
    rownames(out@data) <- rownames(data@data)
    out

  } else {
    out <- data[idx, , drop = FALSE]
    out[, coord.names] <- data[, coord.names]
    if (is.data.frame(out)) rownames(out) <- NULL
    out
  }
}
