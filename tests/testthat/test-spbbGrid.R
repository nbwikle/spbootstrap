# Grid with n.x = 12 columns and n.y = 8 rows (non-square), unit spacing,
# rows deliberately shuffled.
make_grid <- function(n.x = 12, n.y = 8, shuffle = TRUE) {
  df <- expand.grid(x = seq_len(n.x) - 0.5, y = seq_len(n.y) - 0.5)
  df$z <- seq_len(nrow(df))
  if (shuffle) df <- df[sample(nrow(df)), ]
  rownames(df) <- NULL
  df
}

# Check that each b.l x b.l tile of a bootstrap replicate is a rigid,
# un-rotated translation of a contiguous block in the original grid.
expect_block_structure <- function(idx, coords, b.l) {
  src  <- coords[idx, ]
  dx   <- src[, 1] - coords[, 1]
  dy   <- src[, 2] - coords[, 2]
  # all indices valid and every destination cell filled
  expect_false(anyNA(idx))
  expect_true(all(idx >= 1 & idx <= nrow(coords)))
  # neighbouring destination cells in the same tile share the same shift
  n.same <- 0
  for (r in seq_len(nrow(coords))) {
    nb <- which(abs(coords[, 1] - coords[r, 1] - 1) < 1e-8 &
                abs(coords[, 2] - coords[r, 2]) < 1e-8)
    if (length(nb) == 1 && dx[nb] == dx[r] && dy[nb] == dy[r]) n.same <- n.same + 1
  }
  # at most one tile boundary per b.l cells in each row (plus one when the
  # tiling is shifted), so most pairs match
  n.pairs <- sum(sapply(seq_len(nrow(coords)), function(r)
    any(abs(coords[, 1] - coords[r, 1] - 1) < 1e-8 &
        abs(coords[, 2] - coords[r, 2]) < 1e-8)))
  n.rows <- length(unique(coords[, 2]))
  expect_gte(n.same, n.pairs * (1 - 1 / b.l) - n.rows - 1e-8)
}

test_that("data.frame input works on a non-square, shuffled grid", {
  set.seed(1)
  df <- make_grid()
  b  <- spbbGrid(df, coord.names = c("x", "y"), n.boot = 5, block.l = c(2, 3, 4, 5))

  expect_named(b$bootstraps, c("bl.2", "bl.3", "bl.4", "bl.5"))
  for (bl in c(2, 3, 4, 5)) {
    m <- b$bootstraps[[paste0("bl.", bl)]]
    expect_equal(dim(m), c(nrow(df), 5))
    for (k in 1:5) expect_block_structure(m[, k], df[, c("x", "y")], bl)
  }
})

test_that("block.l equal to both grid dimensions reproduces the grid", {
  set.seed(2)
  df <- make_grid(6, 6)
  b  <- spbbGrid(df, c("x", "y"), n.boot = 3, block.l = 6, shift = FALSE)
  expect_equal(b$bootstraps$bl.6, matrix(seq_len(36), 36, 3))
})

test_that("all input types give identical bootstraps", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")

  r <- terra::rast(nrows = 8, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 8,
                   crs = "local")
  terra::values(r) <- seq_len(terra::ncell(r))
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)   # same order as terra cells

  sf.pts <- sf::st_as_sf(df, coords = c("x", "y"))
  sp.pts <- df; sp::coordinates(sp.pts) <- ~ x + y
  sp.pix <- sp.pts; sp::gridded(sp.pix) <- TRUE
  rl     <- raster::raster(r)

  run <- function(d, cn = NULL) {
    set.seed(10)
    spbbGrid(d, coord.names = cn, n.boot = 4, block.l = c(3, 5))$bootstraps
  }
  ref <- run(df, c("x", "y"))

  expect_identical(run(r), ref)
  expect_identical(run(rl), ref)
  expect_identical(run(sf.pts), ref)
  expect_identical(run(sp.pts), ref)
  expect_identical(run(as.matrix(df), c("x", "y")), ref)

  # SpatialPixels may reorder cells, so compare via coordinates instead
  set.seed(10)
  b.pix <- spbbGrid(sp.pix, n.boot = 4, block.l = c(3, 5))
  pix.xy <- b.pix$coords
  for (bl in names(ref)) {
    key.ref <- paste(df$x[ref[[bl]]], df$y[ref[[bl]]])
    key.pix <- paste(pix.xy$x[b.pix$bootstraps[[bl]]], pix.xy$y[b.pix$bootstraps[[bl]]])
    ord <- match(paste(df$x, df$y), paste(pix.xy$x, pix.xy$y))
    expect_identical(matrix(key.pix, ncol = 4)[ord, ], matrix(key.ref, ncol = 4))
  }
})

test_that("raster indices are cell numbers", {
  skip_if_not_installed("terra")
  set.seed(3)
  r <- terra::rast(nrows = 10, ncols = 7, crs = "local")
  terra::values(r) <- rnorm(terra::ncell(r))
  b <- spbbGrid(r, n.boot = 2, block.l = 3)
  expect_equal(nrow(b$bootstraps$bl.3), terra::ncell(r))
  expect_true(all(b$bootstraps$bl.3 %in% seq_len(terra::ncell(r))))
  # resampled layer can be built directly from cell values
  r.boot <- r; terra::values(r.boot) <- terra::values(r)[b$bootstraps$bl.3[, 1]]
  expect_equal(terra::ncell(r.boot), terra::ncell(r))
})

test_that("point process data are accepted as data.frame, sf and sp", {
  skip_if_not_installed("sp")
  df  <- make_grid(10, 6, shuffle = FALSE)
  pts <- data.frame(px = runif(80, 0, 10), py = runif(80, 0, 6))
  pts.sf <- sf::st_as_sf(pts, coords = c("px", "py"))
  pts.sp <- pts; sp::coordinates(pts.sp) <- ~ px + py

  run <- function(p, cn = NULL) {
    set.seed(4)
    spbbGrid(df, c("x", "y"), n.boot = 3, block.l = 2, pt.data = p, pt.coord.names = cn)
  }
  ref <- run(pts, c("px", "py"))
  expect_length(ref$pt.bootstraps$bl.2, 3)
  expect_true(all(unlist(ref$pt.bootstraps) %in% seq_len(nrow(pts))))
  expect_identical(run(pts.sf)$pt.bootstraps, ref$pt.bootstraps)
  expect_identical(run(pts.sp)$pt.bootstraps, ref$pt.bootstraps)
  # the grid bootstrap is unaffected by adding points
  set.seed(4)
  expect_identical(spbbGrid(df, c("x", "y"), n.boot = 3, block.l = 2)$bootstraps,
                   ref$bootstraps)
})

test_that("with evenly dividing blocks, points land in the source blocks' cells", {
  set.seed(5)
  df  <- make_grid(8, 6, shuffle = FALSE)
  pts <- data.frame(px = runif(200, 0, 8), py = runif(200, 0, 6))
  b   <- spbbGrid(df, c("x", "y"), n.boot = 3, block.l = 2, pt.data = pts,
                  pt.coord.names = c("px", "py"))
  cell.of <- function(x, y) paste(floor(x), floor(y))
  pt.cell <- cell.of(pts$px, pts$py)
  for (k in 1:3) {
    src.cells <- cell.of(df$x[b$bootstraps$bl.2[, k]], df$y[b$bootstraps$bl.2[, k]])
    # multiset of resampled points equals points in the resampled source cells
    expected <- unlist(lapply(src.cells, function(cc) which(pt.cell == cc)))
    expect_equal(sort(b$pt.bootstraps$bl.2[[k]]), sort(expected))
  }
})

test_that("invalid inputs give informative errors", {
  df <- make_grid(6, 5)
  expect_error(spbbGrid(df, n.boot = 2, block.l = 2), "coord.names")
  expect_error(spbbGrid(df, c("x", "q"), n.boot = 2, block.l = 2), "not found")
  expect_error(spbbGrid(df[-1, ], c("x", "y"), n.boot = 2, block.l = 2), "complete grid")
  expect_error(spbbGrid(df, c("x", "y"), n.boot = 2, block.l = 6), "between 1 and 5")
  expect_error(spbbGrid(df, c("x", "y"), n.boot = 2), "block.l")
  expect_error(spbbGrid(list(1), n.boot = 2, block.l = 2), "must be a SpatRaster")
  poly <- sf::st_as_sf(sf::st_make_grid(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 3, ymax = 3)), n = 3))
  expect_error(spbbGrid(poly, n.boot = 2, block.l = 2), "POINT geometry")
})


test_that("spbbGrid shift = TRUE moves the tiling by whole cells between replicates", {
  set.seed(17)
  df <- expand.grid(x = 1:12, y = 1:9)
  src <- function(shift) {
    spbbGrid(df, c("x", "y"), n.boot = 40, block.l = 3, shift = shift)$bootstraps$bl.3
  }
  # cells (3, 5) and (4, 5) straddle a fixed tile edge (x = 3 | 4) when shift = FALSE;
  # they are neighbours from the same source block iff their sources are x-neighbours
  a <- which(df$x == 3 & df$y == 5); b <- which(df$x == 4 & df$y == 5)
  together <- function(m) df$x[m[b, ]] - df$x[m[a, ]] == 1 & df$y[m[b, ]] == df$y[m[a, ]]
  expect_false(any(together(src(FALSE))))
  expect_gt(sum(together(src(TRUE))), 10)

  # every cell is still filled, with valid indices
  m <- src(TRUE)
  expect_false(anyNA(m))
  expect_true(all(m >= 1 & m <= nrow(df)))
})
