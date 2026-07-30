trapezoidal_weights <- function(grid) {
  grid <- as.numeric(grid)
  
  if (length(grid) < 2 || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    stop("The functional grid must be finite and strictly increasing.")
  }
  
  number_grid_points <- length(grid)
  grid_differences <- diff(grid)
  weights <- numeric(number_grid_points)
  
  weights[1] <- grid_differences[1] / 2
  weights[number_grid_points] <- grid_differences[number_grid_points - 1] / 2
  
  if (number_grid_points > 2) {
    weights[2:(number_grid_points - 1)] <-
      (grid_differences[-1] +
         grid_differences[-length(grid_differences)]) / 2
  }
  
  weights
}
