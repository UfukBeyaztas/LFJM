polynomial_design <- function(time, degree) {
  if (degree < 1) {
    return(matrix(1, length(time), 1))
  }
  
  sapply(0:degree, function(power) time^power)
}