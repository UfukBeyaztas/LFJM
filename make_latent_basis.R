make_latent_basis <- function(
    grid,
    number_components,
    type = c("B0", "B1", "U")
) {
  type <- match.arg(type)
  weights <- trapezoidal_weights(grid)
  x <- 2 * pi * (grid - min(grid)) / diff(range(grid))
  z <- 2 * (grid - min(grid)) / diff(range(grid)) - 1
  basis <- matrix(0, length(grid), number_components)
  
  for (k in seq_len(number_components)) {
    if (type == "B0") {
      basis[, k] <- if (k %% 2 == 1) {
        sin(ceiling(k / 2) * x)
      } else {
        cos((k / 2) * x)
      }
    } else if (type == "B1") {
      basis[, k] <- if (k %% 2 == 1) {
        cos(ceiling(k / 2) * x)
      } else {
        sin((k / 2 + 1) * x)
      }
    } else {
      basis[, k] <- if (k %% 2 == 1) {
        sin(ceiling(k / 2) * x)
      } else {
        cos((k / 2) * x)
      }
    }
  }
  
  orthogonalize_columns(basis, weights)
}