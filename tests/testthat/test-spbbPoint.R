test_that("spbbPoint uses coord.names rather than columns named x and y", {
  pts <- data.frame(lon = runif(300, 0, 50), lat = runif(300, 0, 50), v = 1:300)
  set.seed(6)
  a <- spbbPoint(pts, c("lon", "lat"), n.boot = 3, block.l = 10,
                   box.x = c(0, 50), box.y = c(0, 50))
  pts.xy <- setNames(pts, c("x", "y", "v"))
  set.seed(6)
  b <- spbbPoint(pts.xy, c("x", "y"), n.boot = 3, block.l = 10,
                   box.x = c(0, 50), box.y = c(0, 50))
  for (k in 1:3) {
    expect_named(a$bl.10[[k]], c("lon", "lat", "v"))
    expect_equal(unname(as.matrix(a$bl.10[[k]])), unname(as.matrix(b$bl.10[[k]])))
  }
})

test_that("spbbPoint keeps points inside a box that does not start at 0", {
  pts <- data.frame(e = runif(400, 100, 140), n = runif(400, -30, 0), id = 1:400)
  set.seed(8)
  b <- spbbPoint(pts, c("e", "n"), n.boot = 4, block.l = 10,
                   box.x = c(100, 140), box.y = c(-30, 0))
  for (k in 1:4) {
    d <- b$bl.10[[k]]
    expect_gt(nrow(d), 0.5 * nrow(pts))
    expect_true(all(d$e >= 100 & d$e <= 140 & d$n >= -30 & d$n <= 0))
    # resampled rows are copies of original rows (only coordinates change)
    expect_true(all(d$id %in% pts$id))
  }
})

test_that("sf and sp inputs give the same bootstraps as a data.frame, in their own class", {
  skip_if_not_installed("sp")
  pts <- data.frame(x = runif(300, 0, 40), y = runif(300, 0, 30), v = 1:300)
  pts.sf  <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 32615)
  pts.sfc <- sf::st_geometry(pts.sf)
  pts.sp  <- pts; sp::coordinates(pts.sp) <- ~ x + y
  pts.spo <- sp::SpatialPoints(pts[, c("x", "y")])

  run <- function(d, cn = NULL) {
    set.seed(12)
    spbbPoint(d, cn, n.boot = 3, block.l = c(5, 10), box.x = c(0, 40), box.y = c(0, 30))
  }
  ref <- run(pts, c("x", "y"))
  a.sf <- run(pts.sf); a.sfc <- run(pts.sfc); a.sp <- run(pts.sp); a.spo <- run(pts.spo)

  for (bl in names(ref)) for (k in 1:3) {
    r <- ref[[bl]][[k]]
    xy <- unname(as.matrix(r[, c("x", "y")]))

    expect_s3_class(a.sf[[bl]][[k]], "sf")
    expect_equal(sf::st_crs(a.sf[[bl]][[k]]), sf::st_crs(pts.sf))
    expect_equal(unname(sf::st_coordinates(a.sf[[bl]][[k]])), xy)
    expect_equal(a.sf[[bl]][[k]]$v, r$v)

    expect_s3_class(a.sfc[[bl]][[k]], "sfc")
    expect_equal(unname(sf::st_coordinates(a.sfc[[bl]][[k]])), xy)

    expect_s4_class(a.sp[[bl]][[k]], "SpatialPointsDataFrame")
    expect_equal(unname(sp::coordinates(a.sp[[bl]][[k]])), xy)
    expect_equal(a.sp[[bl]][[k]]$v, r$v)

    expect_s4_class(a.spo[[bl]][[k]], "SpatialPoints")
    expect_equal(unname(sp::coordinates(a.spo[[bl]][[k]])), xy)
  }
})

test_that("block.l cannot exceed the shorter side of the box", {
  pts <- data.frame(x = runif(50, 0, 40), y = runif(50, 0, 30))
  expect_error(spbbPoint(pts, c("x", "y"), n.boot = 2, block.l = 31,
                         box.x = c(0, 40), box.y = c(0, 30)), "at most 30")
  expect_error(spbbPoint(pts, c("x", "y"), n.boot = 2, block.l = c(5, 35),
                         box.x = c(0, 40), box.y = c(0, 30)), "at most 30")
  expect_error(spbbPoint(pts, c("x", "y"), n.boot = 2, block.l = 0,
                         box.x = c(0, 40), box.y = c(0, 30)), "greater than 0")
  expect_error(spbbPoint(pts, c("x", "y"), n.boot = 2, block.l = 30.01,
                         box.x = c(0, 40), box.y = c(0, 30)), "at most 30")
  # the largest allowed block works
  expect_length(spbbPoint(pts, c("x", "y"), n.boot = 2, block.l = 30,
                          box.x = c(0, 40), box.y = c(0, 30))$bl.30, 2)
})

test_that("other invalid inputs give informative errors or warnings", {
  pts <- data.frame(x = runif(50, 0, 10), y = runif(50, 0, 10))
  expect_error(spbbPoint(pts, n.boot = 2, block.l = 2, box.x = c(0, 10), box.y = c(0, 10)),
               "coord.names")
  expect_error(spbbPoint(pts, c("x", "y"), n.boot = 2, box.x = c(0, 10), box.y = c(0, 10)),
               "block.l")
  expect_error(spbbPoint(pts, c("x", "y"), n.boot = 2, block.l = 2, box.x = 5, box.y = c(0, 10)),
               "box.x")
  expect_warning(spbbPoint(pts, c("x", "y"), n.boot = 1, block.l = 2,
                           box.x = c(0, 5), box.y = c(0, 10)), "outside")
  skip_if_not_installed("terra")
  expect_error(spbbPoint(terra::rast(nrows = 5, ncols = 5), n.boot = 1, block.l = 2),
               "spbbGrid")
})

test_that("real-valued block lengths are used as given", {
  # small box in degree-like units, where integer blocks would be useless
  pts <- data.frame(lon = runif(500, -90.5, -89.7), lat = runif(500, 43.1, 43.7), id = 1:500)
  set.seed(13)
  b <- spbbPoint(pts, c("lon", "lat"), n.boot = 3, block.l = c(0.15, 0.25),
                 box.x = c(-90.5, -89.7), box.y = c(43.1, 43.7))
  expect_named(b, c("bl.0.15", "bl.0.25"))
  for (bl in names(b)) for (d in b[[bl]]) {
    expect_true(all(d$lon >= -90.5 & d$lon <= -89.7 & d$lat >= 43.1 & d$lat <= 43.7))
    expect_gt(nrow(d), 0.5 * nrow(pts))
  }

  # block of length 12.5 on a 25 x 25 box tiles it exactly with 4 blocks, and every
  # resampled point keeps its offset from its source block's corner modulo 12.5
  pts <- data.frame(x = runif(400, 0, 25), y = runif(400, 0, 25), id = 1:400)
  set.seed(14)
  d <- spbbPoint(pts, c("x", "y"), n.boot = 1, block.l = 12.5,
                 box.x = c(0, 25), box.y = c(0, 25), shift = FALSE)$bl.12.5[[1]]
  tile <- paste(d$x %/% 12.5, d$y %/% 12.5)
  expect_setequal(unique(tile), c("0 0", "0 1", "1 0", "1 1"))
  # within each tile the shift from the original location is a single constant
  shift <- cbind(d$x - pts$x[d$id], d$y - pts$y[d$id])
  for (t in unique(tile)) {
    s.t <- shift[tile == t, , drop = FALSE]
    expect_lt(max(apply(s.t, 2, function(v) diff(range(v)))), 1e-10)
  }
})

test_that("the box defaults to the data's bounding box, with a warning", {
  skip_if_not_installed("sp")
  pts <- data.frame(x = runif(300, 510, 530), y = runif(300, 4700, 4712), v = 1:300)
  pts.sf <- sf::st_as_sf(pts, coords = c("x", "y"))
  pts.sp <- pts; sp::coordinates(pts.sp) <- ~ x + y
  bb <- sf::st_bbox(pts.sf)
  expect_equal(unname(sp::bbox(pts.sp)[1, ]), unname(bb[c("xmin", "xmax")]))

  run <- function(d, cn = NULL, ...) {
    set.seed(15)
    spbbPoint(d, cn, n.boot = 2, block.l = 4, ...)
  }
  explicit <- run(pts, c("x", "y"),
                  box.x = unname(bb[c("xmin", "xmax")]), box.y = unname(bb[c("ymin", "ymax")]))

  expect_warning(a <- run(pts, c("x", "y")), "bounding box.*study area is known")
  expect_identical(a, explicit)
  expect_warning(a.sf <- run(pts.sf), "bounding box")
  expect_warning(a.sp <- run(pts.sp), "bounding box")
  for (k in 1:2) {
    expect_equal(unname(sf::st_coordinates(a.sf$bl.4[[k]])),
                 unname(as.matrix(explicit$bl.4[[k]][, c("x", "y")])))
    expect_equal(unname(sp::coordinates(a.sp$bl.4[[k]])),
                 unname(as.matrix(explicit$bl.4[[k]][, c("x", "y")])))
  }

  # a supplied box gives no warning; supplying only one side defaults the other
  expect_no_warning(run(pts, c("x", "y"), box.x = c(510, 530), box.y = c(4700, 4712)))
  expect_warning(run(pts, c("x", "y"), box.x = c(510, 530)), "No `box.y` supplied")
  expect_warning(run(pts, c("x", "y"), box.y = c(4700, 4712)), "No `box.x` supplied")
  expect_warning(run(pts, c("x", "y")), "No `box.x` or `box.y` supplied")
})


test_that("spbbPoint shift = TRUE moves the tiling between replicates", {
  set.seed(16)
  pts <- data.frame(x = runif(3000, 0, 30), y = runif(3000, 0, 30), id = 1:3000)
  across <- function(shift) {
    reps <- spbbPoint(pts, c("x", "y"), n.boot = 20, block.l = 10, box.x = c(0, 30),
                      box.y = c(0, 30), shift = shift, output = "indices")$bl.10
    sapply(reps, function(r) {
      id <- paste(round(r$xy[, 1] - pts$x[r$idx], 8), round(r$xy[, 2] - pts$y[r$idx], 8))
      left  <- r$xy[, 1] > 9.8 & r$xy[, 1] < 10
      right <- r$xy[, 1] > 10 & r$xy[, 1] < 10.2
      length(intersect(id[left], id[right])) > 0
    })
  }
  expect_false(any(across(FALSE)))
  expect_gt(sum(across(TRUE)), 10)
  # shifted replicates still cover the whole box
  r <- spbbPoint(pts, c("x", "y"), n.boot = 1, block.l = 10, box.x = c(0, 30),
                 box.y = c(0, 30))$bl.10[[1]]
  expect_true(all(r$x >= 0 & r$x <= 30 & r$y >= 0 & r$y <= 30))
  expect_gt(nrow(r), 2400)
})
