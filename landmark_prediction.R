landmark_prediction <- function(object, landmark) {
  longitudinal_data <- object$longitudinal_data
  survival_data <- object$survival_data
  subject_ids <- survival_data$ID
  
  retained_history <-
    longitudinal_data$obstime <=
    landmark + 1e-10
  
  history_data <- longitudinal_data[
    retained_history,
    ,
    drop = FALSE
  ]
  history_functional_data <-
    object$functional_data[
      retained_history,
      ,
      drop = FALSE
    ]
  
  if (nrow(history_data) == 0) {
    stop(
      "No longitudinal observations are available at the requested landmark."
    )
  }
  
  history_projection <- predict(
    object$lfpca,
    functional_data =
      history_functional_data,
    subject =
      history_data$ID,
    time =
      history_data$obstime
  )
  
  between_scores_history_long <-
    history_projection$scoreB_subject[
      history_projection$subject_index,
      ,
      drop = FALSE
    ]
  
  fixed_history <-
    object$beta["(Intercept)"] +
    object$beta["obstime"] *
    history_data$obstime +
    object$beta["Q1"] *
    history_data$Q1 +
    as.numeric(
      between_scores_history_long %*%
        object$theta_B
    )
  
  if (length(object$theta_U) > 0) {
    fixed_history <-
      fixed_history +
      as.numeric(
        history_projection$scoreU_visit %*%
          object$theta_U
      )
  }
  
  predicted_random_intercept <-
    history_random_intercept(
      response = history_data$Y,
      fixed_part = fixed_history,
      subject = history_data$ID,
      visit_time = history_data$obstime,
      subject_ids = subject_ids,
      landmark = landmark,
      random_intercept_variance =
        object$sigma_u^2,
      residual_variance =
        object$sigma_e^2
    )
  
  between_scores_landmark <-
    history_projection$scoreB_subject[
      match(
        subject_ids,
        history_projection$ids
      ),
      ,
      drop = FALSE
    ]
  
  current_value <-
    object$beta["(Intercept)"] +
    object$beta["obstime"] *
    landmark +
    object$beta["Q1"] *
    survival_data$Q1 +
    as.numeric(
      between_scores_landmark %*%
        object$theta_B
    ) +
    predicted_random_intercept
  
  risk <-
    object$gamma0 +
    object$gamma_Q *
    survival_data$Q1 +
    object$alpha *
    current_value
  
  data.frame(
    ID = subject_ids,
    landmark = landmark,
    current_value = as.numeric(current_value),
    random_intercept =
      as.numeric(predicted_random_intercept),
    risk = as.numeric(risk),
    stringsAsFactors = FALSE
  )
}