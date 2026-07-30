plot.lfjm <- function(
    x,
    type = c(
      "coefficients",
      "eigenfunctions",
      "eigenvalues"
    ),
    level = c("between", "within"),
    ...
) {
  type <- match.arg(type)
  level <- match.arg(level)
  
  if (type == "eigenfunctions") {
    plot(
      x$lfpca,
      type = "eigenfunctions",
      level = level,
      ...
    )
    return(invisible(x))
  }
  
  if (type == "eigenvalues") {
    plot(
      x$lfpca,
      type = "eigenvalues",
      level = level,
      ...
    )
    return(invisible(x))
  }
  
  coefficient_matrix <- cbind(
    beta_B = x$beta_B,
    beta_U = x$beta_U,
    alpha_beta_B = x$surv_B_current
  )
  
  matplot(
    x$grid,
    coefficient_matrix,
    type = "l",
    lty = c(1, 2, 3),
    lwd = 2,
    col = c("black", "steelblue", "firebrick"),
    xlab = "Functional domain",
    ylab = "Coefficient",
    main = "Estimated functional coefficients",
    ...
  )
  
  legend(
    "topright",
    legend = c(
      expression(hat(beta)[B]),
      expression(hat(beta)[U]),
      expression(hat(alpha) * hat(beta)[B])
    ),
    lty = c(1, 2, 3),
    lwd = 2,
    col = c("black", "steelblue", "firebrick"),
    bty = "n"
  )
  
  invisible(x)
}
