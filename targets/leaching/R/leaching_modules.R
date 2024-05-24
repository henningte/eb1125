#' m70
#'
#' @export
m70 <- function(t, k_1, l_1, alpha_1, s) {

  res <- (1.0 - l_1) / (1.0 + (alpha_1 - 1.0) * k_1 * t)^(1.0 / (alpha_1 - 1.0))
  res[res >= 1.0] <- 1.0 - s
  res

}


#' m71
#'
#' @export
m71 <- function(t, k_1, l_1, s) {

  res <- (1.0 - l_1) * exp(-k_1 * t)
  res[res >= 1.0] <- 1.0 - s
  res

}
