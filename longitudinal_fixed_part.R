longitudinal_fixed_part <- function(object) {
  longitudinal_data <- object$longitudinal_data
  subject_ids <- object$survival_data$ID
  
  between_scores <- object$lfpca$scoreB_subject[
    as.character(subject_ids),
    ,
    drop = FALSE
  ]
  
  between_scores_long <- between_scores[
    match(longitudinal_data$ID, subject_ids),
    ,
    drop = FALSE
  ]
  
  fixed_part <-
    object$beta["(Intercept)"] +
    object$beta["obstime"] *
    longitudinal_data$obstime +
    object$beta["Q1"] *
    longitudinal_data$Q1 +
    as.numeric(
      between_scores_long %*%
        object$theta_B
    )
  
  if (length(object$theta_U) > 0) {
    fixed_part <-
      fixed_part +
      as.numeric(
        object$lfpca$scoreU_visit %*%
          object$theta_U
      )
  }
  
  fixed_part
}