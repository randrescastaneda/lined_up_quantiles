# Growth factor from year-by-year NA gross growth (named numeric vector)
growth_factor <- function(gf, from_year, to_year) {
  if (from_year == to_year) return(1)
  if (to_year > from_year) {

    yrs <- (from_year + 1L):to_year |>
      as.character()

    prod(gf[yrs])
  } else {
    yrs <- (to_year + 1L):from_year |>
      as.character()

    1 / prod(gf[yrs])
  }
}

deflate_vector <- \(x, gf, from_year, to_year) {
  growth_factor(gf = gf,
                from_year = from_year,
                to_year = to_year) * x
}



# Weighted step-quantiles (left-inverse of CDF)
wquantile <- function(x, w, probs) {
  o  <- radixorder(x)
  x  <- x[o]; w <- w[o]
  cw <- fcumsum(w); W <- cw[length(cw)]

  idx <- findInterval(probs * W, c(0, cw)) + 0L
  x[pmax(1L, pmin(idx, length(x)))]
}

# Weighted CDF at a threshold z
fgt0 <- function(x, w, z) fmean(x <= z , w = w)



# Time-proximity weight on the earlier survey
alpha_t <- function(y0, y1, t) {
  (y1 - t) / (y1 - y0)
}




# OT/Wasserstein interpolation (unconstrained)
qinterp_base <- function(x0, w0,
                         x1, w1,
                         g0 = 1, g1 = 1,
                         alpha,
                         ngrid = 1e5L,
                         mu_target = NULL) {
  p   <- (seq_len(ngrid) - 0.5) / ngrid
  Q0  <- fquantile(x0 * g0, w0, p)
  Q1  <- fquantile(x1 * g1, w1, p)
  Qb  <- alpha * Q0 + (1 - alpha) * Q1
  if (!is.null(mu_target)) {
    s <- mu_target / fmean(Qb)
    Qb <- Qb * s
  }
  list(p = p, Q = Qb)
}




