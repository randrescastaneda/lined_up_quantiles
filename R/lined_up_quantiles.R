library(fastverse)
library(ggplot2)
library(patchwork) # for combining plots
source("R/utils.R")
set.seed(123)


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

res <- lapply(years, compare_dists) |>
  rbindlist()

res[]




# Gini plot
res_long_gini <- melt(res, id.vars = "year",
                      measure.vars = c("gini_mixture", "gini_qb", "gini_avg"),
                      variable.name = "method", value.name = "gini")
p_gini <- ggplot(res_long_gini, aes(x = year, y = gini, color = method, linetype = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Gini paths", x = "Year", y = "Gini") +
  theme_minimal()

# Poverty plot
res_long_pov <- melt(res, id.vars = "year",
                     measure.vars = c("H_mix", "H_xqb", "H_target"),
                     variable.name = "method", value.name = "poverty")
p_pov <- ggplot(res_long_pov, aes(x = year, y = poverty, color = method, linetype = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Poverty paths", x = "Year", y = "Poverty") +
  theme_minimal()

# Mean plot
res_long_mean <- melt(res, id.vars = "year",
                      measure.vars = c("mean_mix", "mean_qb", "mean_target"),
                      variable.name = "method", value.name = "mean")
p_mean <- ggplot(res_long_mean, aes(x = year, y = mean, color = method, linetype = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Mean paths", x = "Year", y = "Mean") +
  theme_minimal()

# Combine plots side by side
(p_gini | p_pov | p_mean)


