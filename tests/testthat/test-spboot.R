grid_df <- function(n.x = 12, n.y = 8) {
  df <- expand.grid(x = seq_len(n.x) - 0.5, y = seq_len(n.y) - 0.5)
  df$z <- rnorm(nrow(df)) + df$x / 4
  df$w <- rnorm(nrow(df))
  df
}

test_that("grid replicates keep the original locations and move the values", {
  set.seed(1)
  df <- grid_df()
  seen <- list()
  stat <- function(d) {
    # the statistic always sees the original geometry, in the original order
    expect_identical(d[, c("x", "y")], df[, c("x", "y")])
    seen[[length(seen) + 1]] <<- d$z
    mean(d$z)
  }
  set.seed(2)
  fit <- spboot(df, stat, n.boot = 4, block.l = 3, coord.names = c("x", "y"))

  # the values are those of the resampled rows, placed at the destination cells
  set.seed(2)
  idx <- spbbGrid(df, c("x", "y"), n.boot = 4, block.l = 3)$bootstraps$bl.3
  for (k in 1:4) expect_identical(seen[[k + 1]], df$z[idx[, k]])
  expect_equal(unname(fit$t[, 1]), colMeans(matrix(df$z[idx], ncol = 4)))
  expect_equal(unname(fit$se), sd(fit$t[, 1]))
})

test_that("a block taken from the upper right is relocated to where it is placed", {
  # 2 x 2 grid of 2 x 2 blocks: with block.l = 2 every tile is a whole source block
  df <- expand.grid(x = 1:4, y = 1:4)
  df$src.x <- df$x; df$src.y <- df$y
  reps <- list()
  set.seed(3)
  spboot(df, function(d) { reps[[length(reps) + 1]] <<- d; 0 }, n.boot = 20,
         block.l = 2, coord.names = c("x", "y"), shift = FALSE)
  for (d in reps[-1]) {
    ll <- d[d$x <= 2 & d$y <= 2, ]           # lower-left tile of the replicate
    # the tile holds one contiguous source block, shifted by a single offset
    expect_length(unique(ll$src.x - ll$x), 1)
    expect_length(unique(ll$src.y - ll$y), 1)
  }
  # at least one replicate moved an upper-right block (x, y in 3:4) to the lower left
  moved <- sapply(reps[-1], function(d) {
    ll <- d[d$x <= 2 & d$y <= 2, ]
    all(ll$src.x >= 3 & ll$src.y >= 3)
  })
  expect_true(any(moved))
})

test_that("all grid input types give the same bootstrap distribution", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")
  r <- terra::rast(nrows = 8, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 8,
                   crs = "local")
  set.seed(4)
  terra::values(r) <- rnorm(terra::ncell(r))
  names(r) <- "z"
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  sf.g <- sf::st_as_sf(df, coords = c("x", "y"))
  sp.g <- df; sp::coordinates(sp.g) <- ~ x + y
  rl   <- raster::raster(r)

  run <- function(d, stat, cn = NULL) {
    set.seed(5)
    spboot(d, stat, n.boot = 10, block.l = 3, coord.names = cn)$t
  }
  ref <- run(df, function(d) mean(d$z[1:20]), c("x", "y"))
  expect_equal(run(r,    function(d) mean(terra::values(d)[1:20])), ref)
  expect_equal(run(rl,   function(d) mean(raster::values(d)[1:20])), ref)
  expect_equal(run(sf.g, function(d) mean(d$z[1:20])), ref)
  expect_equal(run(sp.g, function(d) mean(d$z[1:20])), ref)
})

test_that("point replicates are the relocated points from spbbPoint", {
  set.seed(6)
  pts <- data.frame(x = runif(300, 0, 30), y = runif(300, 0, 20), v = rnorm(300))
  stat <- function(d) c(n = nrow(d), mx = mean(d$x))
  set.seed(7)
  fit <- spboot(pts, stat, n.boot = 5, block.l = 5.5, type = "point",
                coord.names = c("x", "y"), box.x = c(0, 30), box.y = c(0, 20))
  set.seed(7)
  reps <- spbbPoint(pts, c("x", "y"), n.boot = 5, block.l = 5.5,
                    box.x = c(0, 30), box.y = c(0, 20))[[1]]
  expect_equal(unname(fit$t), unname(t(sapply(reps, stat))))
  expect_warning(
    spboot(pts, stat, n.boot = 2, block.l = 5, type = "point", coord.names = c("x", "y")),
    "bounding box"
  )
})

test_that("per-element block lengths match separate runs", {
  set.seed(8)
  df <- grid_df()
  stat <- function(d) c(a = mean(d$z), b = mean(d$w))

  set.seed(9)
  both <- spboot(df, stat, n.boot = 6, block.l = c(b = 4, a = 2), coord.names = c("x", "y"))
  set.seed(9)
  a2 <- spboot(df, stat, n.boot = 6, block.l = 2, coord.names = c("x", "y"))
  b4 <- spboot(df, stat, n.boot = 6, block.l = 4, coord.names = c("x", "y"))

  expect_equal(both$block.l, c(a = 2, b = 4))
  expect_equal(both$t[, "a"], a2$t[, "a"])
  expect_equal(both$t[, "b"], b4$t[, "b"])
  expect_equal(both$se, c(a = a2$se[["a"]], b = b4$se[["b"]]))

  expect_error(spboot(df, stat, 2, block.l = c(2, 3, 4), coord.names = c("x", "y")),
               "length 1 or 2")
  expect_error(spboot(df, stat, 2, block.l = c(a = 2, q = 3), coord.names = c("x", "y")),
               "must match")
})

test_that("failed replicates become NA with a warning", {
  set.seed(10)
  df <- grid_df()
  calls <- 0
  flaky <- function(d) {
    calls <<- calls + 1
    if (calls %in% c(3, 5)) stop("did not converge")
    mean(d$z)
  }
  expect_warning(
    fit <- spboot(df, flaky, n.boot = 6, block.l = 2, coord.names = c("x", "y")),
    "failed on 2 replicate.*did not converge"
  )
  expect_equal(sum(is.na(fit$t)), 2)
  expect_equal(fit$n.error[["t1"]], 2L)
  expect_true(is.finite(fit$se))
  expect_error(
    spboot(df, function(d) if (nrow(d) > 0 && identical(d, df)) 1 else c(1, 2),
           n.boot = 2, block.l = 2, coord.names = c("x", "y")),
    "returned 2 value"
  )
})

test_that("apply.fun can evaluate replicates in parallel without changing results", {
  skip_on_os("windows")
  set.seed(11)
  df <- grid_df()
  stat <- function(d) coef(lm(z ~ x + w, d))
  set.seed(12)
  a <- spboot(df, stat, n.boot = 8, block.l = 3, coord.names = c("x", "y"))
  set.seed(12)
  b <- spboot(df, stat, n.boot = 8, block.l = 3, coord.names = c("x", "y"),
              apply.fun = function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 2))
  expect_identical(a$t, b$t)
  expect_named(a$t0, c("(Intercept)", "x", "w"))
})

# ---- blockLength --------------------------------------------------------------

test_that("blockLength applies the plug-in formula to the pilot variances", {
  set.seed(13)
  df <- grid_df(30, 20)
  stat <- function(d) c(m = mean(d$z), s = unname(coef(lm(z ~ w, d))[2]))
  bl <- suppressWarnings(blockLength(df, stat, n.boot = 30, coord.names = c("x", "y")))

  n  <- 600
  b1 <- round(n^(1/4)); b2 <- round(0.5 * n^(1/6))
  expect_equal(unname(bl$pilot), c(b1, b2, 2 * b2))
  expect_equal(bl$var["b1", ], bl$fits$b1$se^2)

  B   <- 2 * b2 * (bl$var["b2", ] - bl$var["b3", ])
  raw <- ceiling((B^2 * n / 2 / bl$var["b1", ]^2)^(1/4) * sqrt(3 / 2))
  expect_equal(bl$block.l, pmin(pmax(raw, 1), 20))
  expect_true(all(bl$block.l == round(bl$block.l)))
  expect_named(bl$block.l, c("m", "s"))

  # result feeds straight into spboot
  fit <- spboot(df, stat, n.boot = 5, block.l = bl$block.l, coord.names = c("x", "y"))
  expect_equal(fit$block.l, bl$block.l)
})

test_that("blockLength for point data works in coordinate units", {
  set.seed(14)
  pts <- data.frame(x = runif(1500, 0, 60), y = runif(1500, 0, 40), v = rnorm(1500))
  stat <- function(d) mean(d$v)
  bl <- suppressWarnings(blockLength(pts, stat, n.boot = 20, type = "point",
                                     coord.names = c("x", "y"),
                                     box.x = c(0, 60), box.y = c(0, 40)))
  unit <- sqrt(60 * 40 / 1500)
  expect_equal(bl$unit, unit)
  expect_equal(unname(bl$pilot), c(1500^(1/4), 0.5 * 1500^(1/6), 1500^(1/6)) * unit)

  # computed in count units, converted once at the end
  b2  <- 0.5 * 1500^(1/6)
  B   <- 2 * b2 * (bl$var["b2", ] - bl$var["b3", ])
  raw <- (B^2 * 1500 / 2 / bl$var["b1", ]^2)^(1/4) * sqrt(3 / 2) * unit
  expect_equal(unname(bl$block.l), unname(pmin(pmax(raw, unit), 40)))

  # one bounding-box warning, not one per pilot run
  w <- character(0)
  withCallingHandlers(
    blockLength(pts, stat, n.boot = 5, type = "point", coord.names = c("x", "y")),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_equal(sum(grepl("bounding box", w)), 1)
})

test_that("blockLength rejects pilot blocks larger than the grid", {
  df <- grid_df(6, 5)
  expect_error(blockLength(df, function(d) mean(d$z), n.boot = 2, coord.names = c("x", "y"),
                           c1 = 3), "reduce `c1` or `c2`")
})

test_that("blockLength keeps names for a single-valued statistic", {
  set.seed(15)
  df <- grid_df(20, 15)
  bl <- suppressWarnings(blockLength(df, function(d) c(avg = mean(d$z)), n.boot = 10,
                                     coord.names = c("x", "y")))
  expect_named(bl$block.l, "avg")
  expect_named(bl$B.hat, "avg")
  fit <- spboot(df, function(d) c(avg = mean(d$z)), n.boot = 3, block.l = bl$block.l,
                coord.names = c("x", "y"))
  expect_named(fit$se, "avg")
})
