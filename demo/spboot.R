### spboot.R
### Demo of blockLength() and spboot(): spatial block bootstrap standard errors
### for (1) a regression slope on a raster, (2) an average treatment effect
### with spatial confounding, and (3) a regression on point data.
###
### Run with: demo("spboot", package = "spbootstrap")
### Takes 2-3 minutes; set `fast <- TRUE` for a quicker, noisier run, or use
### several cores via `apply.fun` (see section 2).

library(spbootstrap)
library(terra)
library(mgcv)
library(sf)

fast   <- FALSE
n.pilot <- if (fast) 30 else 50    # replicates per pilot block length
n.final <- if (fast) 100 else 200  # replicates for the final standard error

set.seed(2026)

# Smooth, spatially correlated random field on the grid of `template`:
# white noise averaged over a width x width moving window, then standardized.
smoothField <- function(template, width) {
  r <- terra::setValues(template, rnorm(terra::ncell(template)))
  r <- terra::focal(r, w = width, fun = "mean", na.rm = TRUE)
  (r - terra::global(r, "mean")[[1]]) / terra::global(r, "sd")[[1]]
}

# 50 x 40 grid of unit cells
grid <- terra::rast(nrows = 40, ncols = 50, xmin = 0, xmax = 50, ymin = 0, ymax = 40,
                    crs = "local")   # plain Cartesian units, not lon/lat


#------------------------------------------------------------------------#
#--- 1. Regression slope on a raster                                  ---#
#------------------------------------------------------------------------#

# Covariate and error are both spatially correlated, so the usual lm() standard
# error, which assumes independent errors, is too small.

simRegression <- function() {
  x1 <- smoothField(grid, width = 9)
  e  <- smoothField(grid, width = 9)
  r  <- c(x1, 1 + 0.5 * x1 + e)
  names(r) <- c("x1", "resp")   # avoid "y", which is also a coordinate name
  r
}

# simulate data
r <- simRegression()
# plot raster data
plot(r)

# The statistic: any function of one dataset that returns a (named) numeric
# vector. Here `d` is a SpatRaster, because the data is a SpatRaster.
slope <- function(d) {
  df <- terra::as.data.frame(d)
  c(slope = unname(coef(lm(resp ~ x1, data = df))["x1"]))
}

# Step 1: choose the block length (in grid cells).
bl.reg <- blockLength(r, slope, n.boot = n.pilot)
bl.reg

# Step 2: bootstrap at the chosen block length.
fit.reg <- spboot(r, slope, n.boot = n.final, block.l = bl.reg$block.l)
fit.reg

# Compare with the naive standard error and with the true sampling standard
# deviation, found by simulating new datasets from the same model.
naive.se <- summary(lm(resp ~ x1, data = terra::as.data.frame(r)))$coefficients["x1", "Std. Error"]
true.se  <- sd(replicate(if (fast) 100 else 500, slope(simRegression())))

round(c(naive.lm = naive.se, block.bootstrap = fit.reg$se[["slope"]], truth = true.se), 4)

# The block bootstrap is far closer to the truth than the naive standard error,
# but still somewhat low. The fields here are correlated over ~9 cells on a
# 50 x 40 grid, so the data hold only a few dozen effectively independent
# pieces; block bootstrap standard errors tend to be too small in that setting.

# What a replicate looks like: the values of blocks drawn from anywhere in the
# grid, placed side by side on the original grid. spbbGrid() returns the cell
# numbers; spboot() does this rearrangement internally.
b <- bl.reg$block.l[["slope"]]
idx <- spbbGrid(r, n.boot = 1, block.l = b)$bootstraps[[1]]
r.boot <- terra::setValues(r, terra::values(r)[idx[, 1], ])

op <- par(mfrow = c(1, 2))
terra::plot(r$resp, main = "Original response", legend = FALSE)
terra::plot(r.boot$resp, main = sprintf("Bootstrap replicate (%d x %d blocks)", b, b),
            legend = FALSE)
par(op)


#------------------------------------------------------------------------#
#--- 2. Average treatment effect with spatial confounding              ---#
#------------------------------------------------------------------------#

# An unmeasured confounder U(x, y), a smooth function of location, raises both
# the chance of treatment and the outcome. Both the propensity score and the
# outcome model adjust for it through s(x, y). (If U were rougher than s(x, y)
# can represent, the estimates would be biased: spatial confounding.)
# Replicates carry their *new* coordinates, so s(x, y) is fit to the bootstrap
# geometry. A smaller 40 x 30 grid keeps the gam fits quick.

grid.ate <- terra::rast(nrows = 30, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 30,
                        crs = "local")

eps <- smoothField(grid.ate, width = 5)
ate.df <- terra::as.data.frame(eps, xy = TRUE)
names(ate.df) <- c("x", "y", "eps")
ate.df$U <- as.numeric(scale(
  sin(ate.df$x / 7) + cos(ate.df$y / 5) + 0.5 * sin((ate.df$x + ate.df$y) / 9)
))
ate.df$treat   <- rbinom(nrow(ate.df), 1, plogis(1.2 * ate.df$U))
ate.df$outcome <- 1 + 2 * ate.df$treat + 3 * ate.df$U + ate.df$eps   # true ATE = 2
ate.df$U <- ate.df$eps <- NULL                                         # U is unobserved

# The statistic returns three estimates; each gets its own block length.
ateEstimates <- function(d, k = 30) {
  ps <- fitted(gam(treat ~ s(x, y, k = k), family = binomial, data = d))
  ps <- pmin(pmax(ps, 0.01), 0.99)
  w1 <- d$treat / ps
  w0 <- (1 - d$treat) / (1 - ps)
  c(
    HT         = mean(w1 * d$outcome) - mean(w0 * d$outcome),
    Hajek      = sum(w1 * d$outcome) / sum(w1) - sum(w0 * d$outcome) / sum(w0),
    regression = unname(coef(gam(outcome ~ treat + s(x, y, k = k), data = d))["treat"])
  )
}

# Extra arguments (here `k`, the basis dimension of s(x, y)) are passed on to the
# statistic. To use several cores, add e.g.
#   apply.fun = function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 4)
bl.ate <- blockLength(ate.df, ateEstimates, n.boot = n.pilot,
                      coord.names = c("x", "y"), k = 30)
bl.ate

# spboot() runs once per distinct block length and keeps the matching estimates.
fit.ate <- spboot(ate.df, ateEstimates, n.boot = n.final, block.l = bl.ate$block.l,
                  coord.names = c("x", "y"), k = 30)
fit.ate

# 95% normal-approximation intervals (true ATE = 2)
cbind(
  estimate = fit.ate$t0,
  lower    = fit.ate$t0 - 1.96 * fit.ate$se,
  upper    = fit.ate$t0 + 1.96 * fit.ate$se
)


#------------------------------------------------------------------------#
#--- 3. Point data (sf)                                                ---#
#------------------------------------------------------------------------#

# 1500 points in a 50 x 40 study area, with a mark depending on a covariate;
# both inherit spatial correlation from smooth fields. Block lengths are now
# in coordinate units and need not be whole numbers.

pts <- data.frame(x = runif(1500, 0, 50), y = runif(1500, 0, 40))
cells <- terra::cellFromXY(grid, as.matrix(pts))
f1 <- smoothField(grid, width = 9)
f2 <- smoothField(grid, width = 9)
pts$elev <- terra::values(f1)[cells] + rnorm(1500, sd = 0.3)
pts$mark <- 2 - 0.8 * pts$elev + terra::values(f2)[cells] + rnorm(1500, sd = 0.3)
pts.sf <- sf::st_as_sf(pts, coords = c("x", "y"))

# Replicates are sf objects with relocated geometry.
elevEffect <- function(d) c(elev = unname(coef(lm(mark ~ elev, data = d))["elev"]))

# The study area is known, so supply it rather than relying on the bounding box.
bl.pt <- blockLength(pts.sf, elevEffect, n.boot = n.pilot, type = "point",
                     box.x = c(0, 50), box.y = c(0, 40))
bl.pt

fit.pt <- spboot(pts.sf, elevEffect, n.boot = n.final, block.l = bl.pt$block.l,
                 type = "point", box.x = c(0, 50), box.y = c(0, 40))
fit.pt

naive.pt <- summary(lm(mark ~ elev, data = pts))$coefficients["elev", "Std. Error"]
round(c(naive.lm = naive.pt, block.bootstrap = fit.pt$se[["elev"]]), 4)

# One replicate: blocks of points placed side by side (dashed lines show tiles).
# Points are coloured by `mark` on a common scale, so the moved patches are visible.
rep1 <- spbbPoint(pts.sf, n.boot = 1, block.l = bl.pt$block.l[["elev"]],
                  box.x = c(0, 50), box.y = c(0, 40))[[1]][[1]]
brks <- seq(min(pts$mark), max(pts$mark), length.out = 21)
markCol <- function(m) hcl.colors(20)[cut(m, brks, include.lowest = TRUE)]

op <- par(mfrow = c(1, 2))
plot(sf::st_geometry(pts.sf), pch = 16, cex = 0.6, col = markCol(pts.sf$mark),
     main = "Original points (colour = mark)", axes = TRUE)
plot(sf::st_geometry(rep1), pch = 16, cex = 0.6, col = markCol(rep1$mark),
     main = sprintf("Bootstrap replicate (%.1f x %.1f blocks)",
                    bl.pt$block.l[["elev"]], bl.pt$block.l[["elev"]]), axes = TRUE)
abline(v = seq(0, 50, by = bl.pt$block.l[["elev"]]),
       h = seq(0, 40, by = bl.pt$block.l[["elev"]]), lty = 2, col = "grey40")
par(op)
