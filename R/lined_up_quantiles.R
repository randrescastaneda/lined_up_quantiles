library(fastverse)


# =================== Simulate surveys & NA path ==================
# means
m0 <- 8
m1 <- 10


n0 <- 8000L; n1 <- 9000L
x0 <- rlnorm(n0, meanlog = 8.9, sdlog = 0.55)   # year y0
# normalize to ne mean
x0 <- (m0/fmean(x0))*x0

x1 <- rlnorm(n1, meanlog = 9.1, sdlog = 0.60)   # year y1
x1 <- (m1/fmean(x1))*x1

w0 <- rexp(n0)
w1 <- rexp(n1)

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

# Rescale the two original surveys to match NA means in their own years
x0 <- x0 * (mu_na[as.character(y0)] / fmean(x0, w0))
x1 <- x1 * (mu_na[as.character(y1)] / fmean(x1, w1))

# Poverty line to track (absolute; PPP etc.)
z <- unname(mu_na[as.character(y0)]) * 0.75

# Time-proximity weight on the earlier survey
alpha_t <- function(t) (y1 - t) / (y1 - y0)
