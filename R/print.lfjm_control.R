print.lfjm_control <- function(x, ...) {
  cat("LFJM control parameters\n")
  cat("  Between-level PVE:", x$pve_between, "\n")
  cat("  Within-level PVE: ", x$pve_within, "\n")
  cat("  Between components:",
      if (is.null(x$components_between)) "selected by PVE" else x$components_between,
      "\n")
  cat("  Within components: ",
      if (is.null(x$components_within)) "selected by PVE" else x$components_within,
      "\n")
  cat("  Mean degree:       ", x$mean_degree, "\n")
  cat("  Quadrature points: ", x$quadrature_points, "\n")
  cat("  Maximum iterations:", x$max_iterations, "\n")
  invisible(x)
}
