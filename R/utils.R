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
fgt0 <- function(x, w = rep(1, length(x)), z) fmean(x <= z , w = w)



# Time-proximity weight on the earlier survey
alpha_t <- function(y0, y1, t) {
  (y1 - t) / (y1 - y0)
}




# OT/ type 7 interpolation (unconstrained)
qinterp_base <- function(x0, w0,
                         x1, w1,
                         alpha,
                         nbins = 1e5L,
                         mu_target = NULL,
                         g0 = 1, g1 = 1
                         ) {
  p   <- (seq_len(nbins) - 0.5) / nbins

  Q0  <- fquantile(x0 * g0, w = w0, probs =  p, names = FALSE)
  Q1  <- fquantile(x1 * g1, w = w1, probs =  p, names = FALSE)

  Qb  <- alpha * Q0 + (1 - alpha) * Q1
  if (!is.null(mu_target)) {
    s <- mu_target / fmean(Qb)
  } else {
    s <- (fmean(x = x0, w = w0)*alpha + fmean(x = x1, w = w1)*(1-alpha))/fmean(Qb)
  }
  Qb <- Qb * s
  list(p = p, Q = Qb)
}




compare_dists <- function(t
                          # x0, x1,
                          # w0, w1,
                          # y0, y1,
                          # gf, z
                          ) {

    x0t <- deflate_vector(x = x0,
                          gf = gf, # it must be named
                          from_year = y0,
                          to_year = t)


    x1t <- deflate_vector(x = x1,
                          gf = gf, # it must be named
                          from_year = y1,
                          to_year = t)


    alpha <- alpha_t(y0 = y0,
                     y1 = y1,
                     t = t)


    # X quantiles interpolated
    xqb <- qinterp_base(x0 = x0t,
                        w0 = w0,
                        x1 = x1t,
                        w1 = w1,
                        alpha = alpha)


    # quantile measures
    H_qb <- fgt0(xqb$Q, z = z)
    mean_qb = fmean(xqb$Q)
    gini_qb = wbpip::md_compute_gini(xqb$Q, rep(1, length(xqb$Q)))

    # mixture vector
    x_mix <- c(x0t,x1t)
    w_mix <- c(w0 * alpha, w1 * (1 - alpha))
    g_mixture <- wbpip::md_compute_gini(x_mix, w_mix)
    H_mix     <- fgt0(x_mix, w_mix, z)
    mean_mix  <- fmean(x_mix, w = w_mix)

    # Target measures
    H0t         <- fgt0(x0t, w0, z)
    H1t         <- fgt0(x1t, w1, z)
    mean_target <- alpha * fmean(x0t, w = w0) + (1 - alpha) * fmean(x1t, w = w1)
    H_target    <- alpha * H0t + (1 - alpha) * H1t
    gini_avg    <- alpha * wbpip::md_compute_gini(x0t, w0) +
      (1 - alpha) * wbpip::md_compute_gini(x1t, w1)



    data.table(
      year = t,
      alpha = alpha,
      gini_mixture = g_mixture,
      gini_qb = gini_qb,
      gini_avg = gini_avg,
      H_mix = H_mix,
      H_xqb = H_qb,
      H_target = H_target,
      mean_mix = mean_mix,
      mean_qb = mean_qb,
      mean_target = mean_target
    )
}
