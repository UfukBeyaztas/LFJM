print.lfjm_simulation <- function(x, ...) {
  cat("Simulated longitudinal functional joint-model data\n")
  cat("  Setting:              ", x$case_type, "\n")
  cat("  Subjects:             ", nrow(x$data.surv), "\n")
  cat("  Longitudinal records: ", nrow(x$data.long), "\n")
  cat("  Functional grid size: ", length(x$grid), "\n")
  cat("  Event rate:           ", sprintf("%.3f", mean(x$data.surv$event)), "\n")
  invisible(x)
}
