print.summary.lfjm <- function(x, digits = 4, ...) {
  cat("Longitudinal Functional Joint Model\n")
  cat("Estimation by LFPCA and adjusted profile h-likelihood\n\n")
  
  cat("Subjects:             ",
      x$number_subjects,
      "\n")
  cat("Longitudinal records: ",
      x$number_longitudinal_records,
      "\n")
  cat("Between components:   ",
      x$components_between,
      "\n")
  cat("Within components:    ",
      x$components_within,
      "\n")
  cat("Adjusted log-likelihood:",
      format(
        x$logLik,
        digits = digits
      ),
      "\n")
  cat("Convergence code:       ",
      x$convergence,
      "\n\n")
  
  cat("Conditional Stage-2 coefficient table:\n")
  printCoefmat(
    as.matrix(x$coefficients),
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  
  invisible(x)
}