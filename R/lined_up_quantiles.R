library(fastverse)
source("R/utils.R")



# =================== Simulate surveys & NA path ==================
# means
m0 <- 8
m1 <- 10


n0 <- 8000L
n1 <- 9000L

x0 <- rlnorm(n0, meanlog = 8.9, sdlog = 0.55)   # year y0
x1 <- rlnorm(n1, meanlog = 9.1, sdlog = 0.60)   # year y1

w0 <- rexp(n0)
w1 <- rexp(n1)

# sort vectors
o0 <- radixorder(x0)
o1 <- radixorder(x1)

x0 <- x0[o0]
x1 <- x1[o1]

w0 <- w0[o0]
w1 <- w1[o1]




# normalize to ne mean
x0 <- (m0/fmean(x0, w = w0))*x0
x1 <- (m1/fmean(x1, w = w1))*x1


y0 <- 2000L
y1 <- 2005L


years <- y0:y1

# Year-by-year national-accounts gross growth factors (example)
gf <- c(1.040, 1.015, 0.985, 1.060, 1.040) |>
  setNames(years[-1])

# Build NA mean path anchored at y0 to survey mean (replace with external NA mean if you have it)
mu_na <- length(years) |>
  numeric() |>
  setNames(years)

mu_na[as.character(y0)] <- fmean(x0, w = w0)
yrs_i <- seq_along((y0 + 1L):y1) + 1
for (yr in yrs_i) {
  mu_na[yr] <- mu_na[yr-1L] * gf[yr-1L]
}


# Poverty line to track (absolute; PPP etc.)
z <- unname(mu_na[as.character(y0)]) * 0.5 # let's use half the mean
z <- 3

# vectors in t


# =================== New vectors ==================

t <- 2002

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



xqb <- qinterp_base(x0 = x0t,
                    w0 = w0,
                    x1 = x1t,
                    w1 = w1,
                    alpha = alpha) # X quantiles interpolated


# compare FGTs
fgt0(xqb$Q, z = z)


# mixture vector
x_mix <- c(x0t,x1t)
w_mix <- c(w0 * alpha, w1 * (1 - alpha))

H0t <- fgt0(x0t, w0, z)
H1t <- fgt0(x1t, w1, z)
H_target <- alpha * H0t + (1 - alpha) * H1t

g_mixture <- wbpip::md_compute_gini(x_mix, w_mix)
H_mix     <- fgt0(x_mix, w_mix, z)
mean_mix  <- fmean(x_mix, w = w_mix)
