print.lfjm <- function(x, ...) {
  cat("Longitudinal Functional Joint Model\n")
  cat("  Estimation:                LFPCA-APHL\n")
  cat("  Subjects:                  ",
      nrow(x$survival_data),
      "\n")
  cat("  Longitudinal records:      ",
      nrow(x$longitudinal_data),
      "\n")
  cat("  Between-level components:  ",
      x$lfpca$KB,
      "\n")
  cat("  Within-level components:   ",
      x$lfpca$KU,
      "\n")
  cat("  APHL:                      ",
      sprintf("%.3f", x$loglik),
      "\n")
  cat("  Optimizer convergence code:",
      x$convergence,
      "\n")
  invisible(x)
}