clamp_values <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}
