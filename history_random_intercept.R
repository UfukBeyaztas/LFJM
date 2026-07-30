history_random_intercept <- function(
    response,
    fixed_part,
    subject,
    visit_time,
    subject_ids,
    landmark,
    random_intercept_variance,
    residual_variance
) {
  retained <- visit_time <= landmark + 1e-10
  subject_factor <- factor(
    subject[retained],
    levels = subject_ids
  )
  
  observation_count <- as.numeric(
    table(subject_factor)
  )
  
  residual_sum <- as.numeric(
    tapply(
      response[retained] -
        fixed_part[retained],
      subject_factor,
      sum
    )
  )
  residual_sum[!is.finite(residual_sum)] <- 0
  
  if (!is.finite(random_intercept_variance) ||
      random_intercept_variance <= 0) {
    random_intercept_variance <- 1e-8
  }
  
  if (!is.finite(residual_variance) ||
      residual_variance <= 0) {
    residual_variance <- 1e-8
  }
  
  predicted_random_intercept <-
    residual_sum /
    (
      observation_count +
        residual_variance /
        random_intercept_variance
    )
  
  predicted_random_intercept[
    observation_count == 0
  ] <- NA_real_
  names(predicted_random_intercept) <-
    subject_ids
  
  predicted_random_intercept
}