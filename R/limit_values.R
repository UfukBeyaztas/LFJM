limit_values <- function(x, limit = 35) {
  pmin(pmax(x, -limit), limit)
}
