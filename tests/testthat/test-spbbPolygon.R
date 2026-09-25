# L-shaped study area with a square hole: concave corners and an interior ring
l_shape <- function(crs = sf::NA_crs_) {
  outer <- rbind(c(0, 0), c(60, 0), c(60, 20), c(25, 20), c(25, 40), c(0, 40), c(0, 0))
  hole  <- rbind(c(5, 5), c(12, 5), c(12, 12), c(5, 12), c(5, 5))
  sf::st_sfc(sf::st_polygon(list(outer, hole)), crs = crs)
}

l_points <- function(n = 1500, crs = sf::NA_crs_) {
  poly <- l_shape(crs)
  xy <- cbind(runif(3 * n, 0, 60), runif(3 * n, 0, 40))
  xy <- xy[inPolygon(xy, poly), ][1:n, ]
  data.frame(x = xy[, 1], y = xy[, 2], v = rnorm(n), id = seq_len(n))
}

test_that("the valid-corner region matches an exact containment test", {
  set.seed(1)
  poly <- l_shape()
  for (b in c(3, 8, 15)) {
    v  <- validCorners(poly, b)
    cx <- runif(5000, -2, 60); cy <- runif(5000, -2, 40)
    exact <- seq_along(cx) %in% sf::st_contains(poly, makeSquares(cx, cy, b, sf::NA_crs_))[[1]]
    expect_identical(inPolygon(cbind(cx, cy), v), exact)
  }
  expect_error(sampleSquares(poly, 25, 5), "No square of side 25")
})

test_that("source corners are uniform over the valid region", {
  set.seed(2)
  poly <- l_shape()
  xy <- sampleSquares(poly, 5, 20000)
  v  <- validCorners(poly, 5)
  expect_true(all(inPolygon(xy, v)))
  # share of draws in the lower arm should match its share of the valid area
  lower <- sf::st_intersection(v, sf::st_as_sfc(sf::st_bbox(c(xmin = -1, ymin = -1, xmax = 61, ymax = 15))))
  p <- as.numeric(sf::st_area(lower) / sf::st_area(v))
  expect_equal(mean(xy[, 2] < 15), p, tolerance = 0.02)
})

test_that("replicates are rigid block translations that stay inside the boundary", {
  set.seed(3)
  d <- l_points()
  poly <- l_shape()
  b <- 6
  reps <- polygonReplicates(d, poly, n.boot = 5, block.l = b, coord.names = c("x", "y"),
                            shift = FALSE)[[1]]
  for (r in reps) {
    expect_true(all(inPolygon(r$xy, poly)))
    # tile of each moved point, and its shift from where it came from
    tile  <- paste(floor(r$xy[, 1] / b), floor(r$xy[, 2] / b))
    shift <- r$xy - as.matrix(d[r$idx, c("x", "y")])
    for (t in unique(tile)) {
      s.t <- shift[tile == t, , drop = FALSE]
      expect_lt(max(apply(s.t, 2, function(z) diff(range(z)))), 1e-9)
      # the source square (tile corner minus shift) lies inside the boundary
      corner <- as.numeric(strsplit(t, " ")[[1]]) * b - s.t[1, ]
      expect_true(inPolygon(matrix(corner, 1), validCorners(poly, b)))
    }
  }
})

test_that("output keeps the class, attributes and CRS of the input", {
  skip_if_not_installed("sp")
  set.seed(4)
  d <- l_points(600)
  poly <- l_shape(32615)
  d.sf <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
  d.sp <- sf::as_Spatial(d.sf)

  run <- function(x, cn = NULL) {
    set.seed(5)
    spbbPolygon(x, poly, n.boot = 2, block.l = c(5, 8.5), coord.names = cn)
  }
  a.df <- run(d, c("x", "y")); a.sf <- run(d.sf); a.sp <- run(d.sp)
  expect_named(a.df, c("bl.5", "bl.8.5"))
  for (bl in names(a.df)) for (k in 1:2) {
    xy <- unname(as.matrix(a.df[[bl]][[k]][, c("x", "y")]))
    expect_s3_class(a.sf[[bl]][[k]], "sf")
    expect_equal(sf::st_crs(a.sf[[bl]][[k]]), sf::st_crs(32615))
    expect_equal(unname(sf::st_coordinates(a.sf[[bl]][[k]])), xy)
    expect_equal(a.sf[[bl]][[k]]$id, a.df[[bl]][[k]]$id)
    expect_s4_class(a.sp[[bl]][[k]], "SpatialPointsDataFrame")
    expect_equal(unname(sp::coordinates(a.sp[[bl]][[k]])), xy)
  }
})

test_that("boundary handling: default hull, CRS transform, unions and errors", {
  set.seed(6)
  d <- l_points(400)
  d.sf <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)

  expect_warning(spbbPolygon(d, n.boot = 1, block.l = 4, coord.names = c("x", "y")),
                 "convex hull")

  # a boundary in another CRS is transformed to the data's CRS
  poly <- l_shape(32615)
  set.seed(7); a <- spbbPolygon(d.sf, poly, n.boot = 1, block.l = 5)
  set.seed(7); b <- spbbPolygon(d.sf, sf::st_transform(poly, 3857), n.boot = 1, block.l = 5)
  expect_equal(nrow(a[[1]][[1]]), nrow(b[[1]][[1]]), tolerance = 0.05)

  # several features are unioned into one study area
  halves <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(30, 0), c(30, 40), c(0, 40), c(0, 0)))),
    sf::st_polygon(list(rbind(c(30, 0), c(60, 0), c(60, 40), c(30, 40), c(30, 0))))
  ))
  r <- spbbPolygon(d, halves, n.boot = 1, block.l = 20, coord.names = c("x", "y"))[[1]][[1]]
  expect_gt(nrow(r), 0)

  expect_error(spbbPolygon(d, l_shape(), n.boot = 1, block.l = 45, coord.names = c("x", "y")),
               "at most 40")
  expect_error(spbbPolygon(d, l_shape(), n.boot = 1, block.l = 30, coord.names = c("x", "y")),
               "No square of side 30")
  expect_error(spbbPolygon(d, sf::st_sfc(sf::st_point(c(1, 1))), n.boot = 1, block.l = 2,
                           coord.names = c("x", "y")), "POLYGON")
  expect_warning(
    spbbPolygon(rbind(d, data.frame(x = 40, y = 30, v = 0, id = 0)), l_shape(),
                n.boot = 1, block.l = 4, coord.names = c("x", "y")),
    "1 of 401 points lie outside"
  )
})

test_that("spboot and blockLength support type = 'polygon'", {
  set.seed(8)
  d <- l_points(1200)
  poly <- l_shape()
  stat <- function(x) c(mean = mean(x$v), n = nrow(x))

  set.seed(9)
  fit <- spboot(d, stat, n.boot = 4, block.l = 6, type = "polygon",
                coord.names = c("x", "y"), boundary = poly)
  set.seed(9)
  reps <- spbbPolygon(d, poly, n.boot = 4, block.l = 6, coord.names = c("x", "y"))[[1]]
  expect_equal(unname(fit$t), unname(t(sapply(reps, stat))))
  expect_equal(fit$type, "polygon")

  bl <- suppressWarnings(blockLength(d, function(x) mean(x$v), n.boot = 5, type = "polygon",
                                     coord.names = c("x", "y"), boundary = poly))
  area <- as.numeric(sf::st_area(poly))
  expect_equal(bl$unit, sqrt(area / 1200))
  expect_equal(unname(bl$pilot), c(1200^(1/4), 0.5 * 1200^(1/6), 1200^(1/6)) * bl$unit)

  # the convex-hull warning is given once, not once per pilot run
  w <- character(0)
  withCallingHandlers(
    blockLength(d, function(x) mean(x$v), n.boot = 2, type = "polygon", coord.names = c("x", "y")),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_equal(sum(grepl("convex hull", w)), 1)
})

test_that("lon/lat data work in planar degrees, without s2 errors", {
  skip_if_not(sf::sf_use_s2())
  set.seed(10)
  # study area in UTM (as for a state outline), data in lon/lat
  poly.ll <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(-96, 40.5), c(-91, 40.5), c(-91, 42.2), c(-93, 42.2), c(-93, 43.5),
    c(-96, 43.5), c(-96, 43.5), c(-96, 40.5)          # repeated vertex, as in real outlines
  ))), crs = 4326)
  poly.utm <- sf::st_transform(poly.ll, 26915)
  xy <- cbind(runif(3000, -96, -91), runif(3000, 40.5, 43.5))
  xy <- xy[inPolygon(xy, sf::st_set_crs(poly.ll, NA)), ][1:800, ]
  d <- sf::st_as_sf(data.frame(lon = xy[, 1], lat = xy[, 2], v = rnorm(800)),
                    coords = c("lon", "lat"), crs = 4326)

  reps <- suppressWarnings(spbbPolygon(d, poly.utm, n.boot = 3, block.l = 1))
  expect_length(reps$bl.1, 3)
  for (r in reps$bl.1) {
    expect_equal(sf::st_crs(r), sf::st_crs(4326))
    expect_gt(nrow(r), 0)
  }

  # the unit length is in degrees, from the planar area in square degrees
  bl <- suppressWarnings(blockLength(d, function(x) mean(x$v), n.boot = 3,
                                     type = "polygon", boundary = poly.utm))
  planar.area <- as.numeric(sf::st_area(resolveBoundary(poly.utm, d, xy)))
  expect_equal(bl$unit, sqrt(planar.area / 800))
  expect_lt(bl$unit, 1)
})

test_that("a block length in the wrong units gives a clear error, not a memory blow-up", {
  d <- l_points(200)
  big <- data.frame(x = d$x * 1e4, y = d$y * 1e4, v = d$v)      # "metres"
  poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(6e5, 0), c(6e5, 4e5), c(0, 4e5), c(0, 0)))))
  expect_error(spbbPolygon(big, poly, n.boot = 1, block.l = 1, coord.names = c("x", "y")),
               "units of the coordinates")
  expect_error(spbbPoint(big, c("x", "y"), n.boot = 1, block.l = 1,
                         box.x = c(0, 6e5), box.y = c(0, 4e5)), "units of the coordinates")
})


test_that("shift = TRUE moves the tiling between replicates; shift = FALSE fixes it", {
  set.seed(11)
  d <- l_points(4000)
  poly <- l_shape()
  b <- 6

  # For points placed just either side of x = 6 (a tile edge when shift = FALSE),
  # count replicates in which some such pair came from the same source block.
  same.block.across.6 <- function(shift) {
    reps <- polygonReplicates(d, poly, n.boot = 20, block.l = b,
                              coord.names = c("x", "y"), shift = shift)[[1]]
    sapply(reps, function(r) {
      shift.id <- paste(round(r$xy[, 1] - d$x[r$idx], 8), round(r$xy[, 2] - d$y[r$idx], 8))
      left  <- r$xy[, 1] > 5.8 & r$xy[, 1] < 6 & r$xy[, 2] > 20 & r$xy[, 2] < 40
      right <- r$xy[, 1] > 6 & r$xy[, 1] < 6.2 & r$xy[, 2] > 20 & r$xy[, 2] < 40
      length(intersect(shift.id[left], shift.id[right])) > 0
    })
  }
  expect_false(any(same.block.across.6(FALSE)))
  expect_gt(sum(same.block.across.6(TRUE)), 10)

  # with shift = TRUE each tile is still one rigid block inside the boundary
  reps <- polygonReplicates(d, poly, n.boot = 5, block.l = b, coord.names = c("x", "y"))[[1]]
  for (r in reps) {
    expect_true(all(inPolygon(r$xy, poly)))
    delta <- r$xy - as.matrix(d[r$idx, c("x", "y")])
    grp <- paste(round(delta[, 1], 8), round(delta[, 2], 8))
    span <- sapply(split(seq_along(grp), grp), function(i)
      max(diff(range(r$xy[i, 1])), diff(range(r$xy[i, 2]))))
    expect_true(all(span < b))
  }
})
