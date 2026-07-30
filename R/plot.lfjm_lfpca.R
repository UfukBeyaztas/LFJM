plot.lfjm_lfpca <- function(
    x,
    type = c("eigenfunctions", "eigenvalues"),
    level = c("between", "within"),
    ...
) {
  type <- match.arg(type)
  level <- match.arg(level)
  
  if (type == "eigenvalues") {
    eigenvalues <- if (level == "between") {
      x$eigB
    } else {
      x$eigU
    }
    
    if (length(eigenvalues) == 0) {
      plot.new()
      title(main = paste("No", level, "components retained"))
      return(invisible(x))
    }
    
    plot(
      seq_along(eigenvalues),
      eigenvalues,
      type = "b",
      xlab = "Component",
      ylab = "Eigenvalue",
      main = paste(
        toTitleCase(level),
        "level eigenvalues"
      ),
      ...
    )
  } else {
    eigenfunctions <- if (level == "between") {
      x$phi0
    } else {
      x$phiU
    }
    
    if (ncol(eigenfunctions) == 0) {
      plot.new()
      title(main = paste("No", level, "components retained"))
      return(invisible(x))
    }
    
    matplot(
      x$grid,
      eigenfunctions,
      type = "l",
      lty = 1,
      xlab = "Functional domain",
      ylab = "Eigenfunction",
      main = paste(
        toTitleCase(level),
        "level eigenfunctions"
      ),
      ...
    )
  }
  
  invisible(x)
}
