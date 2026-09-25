# ---- warnings for unprojected longitude/latitude input -------------------------

lonlat_pts <- function(n = 400) {
  data.frame(lon = runif(n, -90, -85), lat = runif(n, 40, 44), v = rnorm(n))
}

count_lonlat_warnings <- function(expr) {
  n <- 0
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("longitude/latitude", conditionMessage(w))) n <<- n + 1
    invokeRestart("muffleWarning")
  })
  n
}

test_that("isLonLat recognises geographic CRSs for every input type", {
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  set.seed(1)
  d <- lonlat_pts(20)
  p.ll   <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
  p.proj <- sf::st_transform(p.ll, 5070)
  p.none <- sf::st_as_sf(d, coords = c("lon", "lat"))
  expect_true(isLonLat(p.ll))
  expect_true(isLonLat(sf::st_geometry(p.ll)))
  expect_false(isLonLat(p.proj))
  expect_false(isLonLat(p.none))          # no CRS: cannot tell, no warning
  expect_false(isLonLat(d))               # data.frame: CRS unknown

  s.ll <- sf::as_Spatial(p.ll)
  expect_true(isLonLat(s.ll))
  expect_false(isLonLat(sf::as_Spatial(p.proj)))

  r.ll    <- terra::rast(nrows = 4, ncols = 5, xmin = -90, xmax = -85, ymin = 40, ymax = 44,
                         crs = "EPSG:4326")
  r.local <- terra::rast(nrows = 4, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 4,
                         crs = "local")
  expect_true(isLonLat(r.ll))
  expect_false(isLonLat(r.local))
  expect_true(isLonLat(raster::raster(r.ll)))
})

test_that("each exported function warns once for lon/lat data, and not for projected data", {
  set.seed(2)
  d  <- lonlat_pts()
  ll <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
  pr <- sf::st_transform(ll, 5070)
  stat <- function(x) mean(x$v)

  expect_warning(
    spbbPoint(ll, n.boot = 1, block.l = 1, box.x = c(-90, -85), box.y = c(40, 44)),
    "equal-area or\\s+local coordinate system"
  )
  expect_equal(count_lonlat_warnings(
    spboot(ll, stat, n.boot = 3, block.l = 1, type = "point",
           box.x = c(-90, -85), box.y = c(40, 44))), 1)
  expect_equal(count_lonlat_warnings(
    blockLength(ll, stat, n.boot = 3, type = "point",
                box.x = c(-90, -85), box.y = c(40, 44))), 1)

  bb <- sf::st_bbox(pr)
  expect_equal(count_lonlat_warnings(
    blockLength(pr, stat, n.boot = 3, type = "point",
                box.x = bb[c("xmin", "xmax")], box.y = bb[c("ymin", "ymax")])), 0)

  # the flag is reset after an error, so the next call still checks
  expect_error(suppressWarnings(spbbPoint(ll, n.boot = 1, block.l = 100,
                                          box.x = c(-90, -85), box.y = c(40, 44))))
  expect_equal(count_lonlat_warnings(
    spbbPoint(ll, n.boot = 1, block.l = 1, box.x = c(-90, -85), box.y = c(40, 44))), 1)
})

test_that("gridded lon/lat data and lon/lat point data in spbbGrid both warn", {
  skip_if_not_installed("terra")
  set.seed(3)
  r <- terra::rast(nrows = 8, ncols = 10, xmin = -90, xmax = -85, ymin = 40, ymax = 44,
                   crs = "EPSG:4326")
  terra::values(r) <- rnorm(terra::ncell(r))
  expect_equal(count_lonlat_warnings(spbbGrid(r, n.boot = 1, block.l = 2)), 1)
  expect_equal(count_lonlat_warnings(
    spboot(r, function(x) terra::global(x, "mean")[[1]], n.boot = 2, block.l = 2)), 1)

  pts <- sf::st_as_sf(lonlat_pts(30), coords = c("lon", "lat"), crs = 4326)
  msgs <- character(0)
  withCallingHandlers(spbbGrid(r, n.boot = 1, block.l = 2, pt.data = pts),
                      warning = function(w) { msgs <<- c(msgs, conditionMessage(w))
                                              invokeRestart("muffleWarning") })
  expect_true(any(grepl("^`data`", msgs)))
  expect_true(any(grepl("^`pt.data`", msgs)))
})
