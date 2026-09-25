test_that("output = 'indices' + buildReplicate reproduces output = 'objects' for points", {
  skip_if_not_installed("sp")
  set.seed(1)
  pts <- data.frame(x = runif(500, 0, 40), y = runif(500, 0, 30), v = rnorm(500), id = 1:500)
  pts.sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 32615)
  pts.sp <- sf::as_Spatial(pts.sf)
  poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(40, 0), c(40, 30), c(10, 30), c(0, 0)))),
                     crs = 32615)

  for (d in list(list(pts, c("x", "y")), list(pts.sf, NULL), list(pts.sp, NULL))) {
    x <- d[[1]]; cn <- d[[2]]

    set.seed(2); obj <- spbbPoint(x, cn, n.boot = 3, block.l = c(5, 7.5),
                                  box.x = c(0, 40), box.y = c(0, 30))
    set.seed(2); ind <- spbbPoint(x, cn, n.boot = 3, block.l = c(5, 7.5),
                                  box.x = c(0, 40), box.y = c(0, 30), output = "indices")
    expect_named(ind, names(obj))
    for (bl in names(obj)) for (k in 1:3) {
      expect_named(ind[[bl]][[k]], c("idx", "xy"))
      expect_equal(buildReplicate(x, ind[[bl]][[k]], cn), obj[[bl]][[k]])
    }

    set.seed(3); obj <- suppressWarnings(spbbPolygon(x, poly, n.boot = 2, block.l = 6, coord.names = cn))
    set.seed(3); ind <- suppressWarnings(spbbPolygon(x, poly, n.boot = 2, block.l = 6, coord.names = cn,
                                                     output = "indices"))
    for (k in 1:2) expect_equal(buildReplicate(x, ind$bl.6[[k]], cn), obj$bl.6[[k]])
  }
})

test_that("buildReplicate turns spbbGrid indices into relocated objects of every class", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")
  r <- terra::rast(nrows = 6, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 6, crs = "local")
  set.seed(4)
  terra::values(r) <- rnorm(terra::ncell(r))
  names(r) <- "z"
  df   <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  d.sf <- sf::st_as_sf(df, coords = c("x", "y"))
  d.sp <- df; sp::coordinates(d.sp) <- ~ x + y

  set.seed(5)
  idx <- spbbGrid(df, c("x", "y"), n.boot = 2, block.l = 3)$bootstraps$bl.3[, 2]
  moved.z <- df$z[idx]

  # values move, locations stay: cell i of the replicate holds the value of cell idx[i]
  b.df <- buildReplicate(df, idx, c("x", "y"))
  expect_equal(b.df[, c("x", "y")], df[, c("x", "y")])
  expect_equal(b.df$z, moved.z)

  b.sf <- buildReplicate(d.sf, idx)
  expect_s3_class(b.sf, "sf")
  expect_equal(sf::st_coordinates(b.sf), sf::st_coordinates(d.sf))
  expect_equal(b.sf$z, moved.z)

  b.sp <- buildReplicate(d.sp, idx)
  expect_s4_class(b.sp, "SpatialPointsDataFrame")
  expect_equal(b.sp$z, moved.z)

  b.r <- buildReplicate(r, idx)
  expect_s4_class(b.r, "SpatRaster")
  expect_equal(terra::values(b.r)[, 1], moved.z)

  b.rl <- buildReplicate(raster::raster(r), idx)
  expect_s4_class(b.rl, "RasterLayer")
  expect_equal(raster::values(b.rl), moved.z)
})

test_that("compact replicates survive saveRDS and are much smaller", {
  set.seed(6)
  pts <- data.frame(x = runif(2000, 0, 50), y = runif(2000, 0, 50),
                    matrix(rnorm(2000 * 20), ncol = 20))
  pts.sf <- sf::st_as_sf(pts, coords = c("x", "y"))
  set.seed(7)
  ind <- spbbPoint(pts.sf, n.boot = 5, block.l = 10, box.x = c(0, 50), box.y = c(0, 50),
                   output = "indices")$bl.10
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  saveRDS(ind, f)
  k3 <- buildReplicate(pts.sf, readRDS(f)[[3]])

  set.seed(7)
  obj <- spbbPoint(pts.sf, n.boot = 5, block.l = 10, box.x = c(0, 50), box.y = c(0, 50))$bl.10
  expect_equal(k3, obj[[3]])
  expect_lt(as.numeric(object.size(ind)), as.numeric(object.size(obj)) / 5)
})

test_that("buildReplicate rejects replicates that do not fit the data", {
  df <- expand.grid(x = 1:4, y = 1:3); df$z <- 1:12
  pts <- data.frame(x = 1:5, y = 1:5)
  expect_error(buildReplicate(df, 1:10, c("x", "y")), "one index per cell")
  expect_error(buildReplicate(df, c(1:11, 13L), c("x", "y")), "between 1 and 12")
  expect_error(buildReplicate(df, 1:12), "coord.names")
  expect_error(buildReplicate(pts, list(idx = 1:2), c("x", "y")), "`idx` and `xy`")
  expect_error(buildReplicate(pts, list(idx = 1:2, xy = matrix(0, 3, 2)), c("x", "y")),
               "one row per element")
  expect_error(buildReplicate(pts, list(idx = c(1, 9), xy = matrix(0, 2, 2)), c("x", "y")),
               "between 1 and 5")
  expect_error(buildReplicate(pts, "a", c("x", "y")), "must be list")
})
