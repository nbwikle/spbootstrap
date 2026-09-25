# spbootstrap

Spatial block bootstrap samples for R.

| Function | Data | Returns |
|---|---|---|
| `spbbGrid()` | regular grid (optionally plus off-grid points) | row indices of each bootstrap replicate, per block length |
| `spbbPoint()` | point data (data.frame, sf, sp) in a rectangular box | bootstrapped copies of `data` (same class) with shifted coordinates |
| `spsubGrid()` | regular grid (optionally plus off-grid points) | row indices of overlapping subsampling windows |
| `spbbPolygon()` | point data (data.frame, sf, sp) in an irregular study area (polygon) | bootstrapped copies of `data` (same class) with shifted coordinates |
| `spbb()` | irregular `sf` points | bootstrapped `sf` objects (original implementation; superseded by `spbbPolygon()`) |
| `buildReplicate()` | a compact replicate from any sampler | the full replicate, in the class of the data |
| `spboot()` | any of the above, plus a `statistic` function | replicate estimates and bootstrap standard errors |
| `blockLength()` | same as `spboot()` | Nordman–Lahiri plug-in block length for each element of the statistic |

## Bootstrapping an estimator

Write a function that computes your estimate from one dataset; `spboot()` calls
it on every replicate, with resampled observations relocated to their new
positions (so spatial terms such as `s(x, y)` use the bootstrap geometry).

```r
slope <- function(d) coef(lm(y ~ elevation + precip, data = d))["elevation"]

bl  <- blockLength(df, slope, n.boot = 100, coord.names = c("x", "y"))
fit <- spboot(df, slope, n.boot = 500, block.l = bl$block.l, coord.names = c("x", "y"))
fit$se
```

A statistic can return several values (e.g. `c(HT = ..., Hajek = ..., aIPTW = ...)`);
`blockLength()` then chooses a block length for each, and `spboot()` runs once
per distinct block length.

```r
# install.packages("devtools")
devtools::install_github("nbwikle/spbootstrap")
```

For point data in an irregular study area (e.g. monitoring sites within a state),
use `type = "polygon"` and pass the study area as `boundary`:

```r
fit <- spboot(wells, my_stat, n.boot = 500, block.l = 30000,
              type = "polygon", boundary = iowa)
```

Project longitude/latitude data first (e.g. `sf::st_transform(x, 26915)`), so
block lengths are in metres.

## Replicates for cluster jobs

`spbbPoint()` and `spbbPolygon()` can return replicates in compact form
(`output = "indices"`: the resampled rows and their new coordinates), and
`spbbGrid()` always returns row/cell indices. Save these once, then build one
full replicate per job with `buildReplicate()`, which returns an object of the
same class as the data (sf, sp, terra, raster or data.frame):

```r
reps <- spbbPolygon(wells, iowa, n.boot = 500, block.l = 30000, output = "indices")
saveRDS(reps$bl.30000, "reps.rds")

## in cluster job k
wells.k <- buildReplicate(wells, readRDS("reps.rds")[[k]])
est.k   <- my_stat(wells.k)
```

For grids, pass one column of the index matrix:
`buildReplicate(r, spbbGrid(r, n.boot = 500, block.l = 5)$bootstraps$bl.5[, k])`.

## Tile placement

By default (`shift = TRUE`) the tiling that resampled blocks are placed into is
moved by a random offset for each replicate, so tile edges do not fall in the
same places every time. Use `shift = FALSE` for the classical fixed tiling
(tiles start at the corner of the grid, box or boundary's bounding box), e.g. to
reproduce results from earlier versions.
