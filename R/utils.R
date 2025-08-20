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
