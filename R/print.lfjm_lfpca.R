print.lfjm_lfpca <- function(x, ...) {
  cat("Longitudinal functional principal component analysis\n")
  cat("  Subjects:                 ", length(x$ids), "\n")
  cat("  Functional grid size:     ", length(x$grid), "\n")
  cat("  Between-level components: ", x$KB, "\n")
  cat("  Within-level components:  ", x$KU, "\n")
  cat("  Curve presmoothing:       ",
      if (x$presmooth) "yes" else "no",
      "\n")
  invisible(x)
}
