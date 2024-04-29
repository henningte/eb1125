#' Uncertainty analysis following Eriksson et al. (2019)
#'
#' @param y Variable for which to perform an uncertainty analysis.
#'
#' @param x Variable conditional on which the uncertainty analysis is performed.
#'
#' @param nbin Positive integer value. Number of bins into which `x` is cut.
#'
#' @export
leaching_gsa_eriksson2019 <- function(y, x, nbin = 10) {

  # normalize
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x - x_min, na.rm = TRUE)
  x_norm <- (x - x_min) / x_max

  # define bins
  bins <- seq(from = 0, to = 1, length.out = nbin + 1L)
  bins <-
    tibble::tibble(
      lower = bins[-length(bins)],
      upper = bins[-1],
      mean = (lower + upper) / 2
    )

  # template with unconditional values
  res <-
    tibble::tibble(
      id_sample = seq_along(x),
      y_mean = mean(y, na.rm = TRUE)
    )

  # loop over bins and compute output
  res <-
    purrr::map_dfr(seq_len(nrow(bins)), function(i) {
      index <-
        if(i == nbin) {
          x_norm >= bins$lower[[i]] & x_norm <= bins$upper[[i]]
        } else {
          x_norm >= bins$lower[[i]] & x_norm < bins$upper[[i]]
        }
      bins %>%
        dplyr::slice(i) %>%
        dplyr::mutate(
          res_bin =
            list(
              res %>%
                dplyr::mutate(
                  y_mean_bin =
                    purrr::map_dbl(seq_along(y), function(j) {
                      if(sum(index[[j]]) == 0) {
                        NA_real_
                      } else {
                        mean(y[[j]][index[[j]]], na.rm = TRUE)
                      }
                    }),
                  y = y,
                  r_j = sum(index)
                )
            )
        )

    }) %>%
    tidyr::unnest(res_bin)


  res %>%
    dplyr::group_by(id_sample) %>%
    dplyr::summarise(
      y_var_bin = 1 / sum(r_j) * sum(r_j * (y_mean_bin - y_mean)^2, na.rm = TRUE),
      y_var = 1 / sum(r_j) * sum((y[[1]] - y_mean[[1]])^2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      s_i = y_var_bin / y_var
    )

}



#' Function to perform the uncertainty analysis
#'
#' @param x_stan_draws .
#'
#' @param nbin .
#'
#' @param data_template A data frame to which the computed sensitivity analysis
#' results are bound. Must have at least the column `id_sample`.
#'
#' @param ... .
#'
#' @export
leaching_make_uncertainty_analysis_1 <- function(x_stan_draws, nbin = 10, data_template, ...) {

  .dots <- list(...)

  # compute sensitivity indices
  res <-
    x_stan_draws %>%
    dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
    dplyr::filter(! is.na(k_2)) %>%
    dplyr::mutate(
      s_k_2_l_2 =
        if(.dots$has_leaching) {
          leaching_gsa_eriksson2019(y = k_2, x = l_2, nbin = nbin) %>%
            dplyr::pull(s_i)
        } else {
          NA_real_
        },
      s_k_2_alpha_2 =
        if(.dots$has_alpha) {
          leaching_gsa_eriksson2019(y = k_2, x = alpha_2, nbin = nbin) %>%
            dplyr::pull(s_i)
        } else {
          NA_real_
        },
      s_l_2_alpha_2 =
        if(.dots$has_alpha & .dots$has_leaching) {
          leaching_gsa_eriksson2019(y = k_2, x = alpha_2, nbin = nbin) %>%
            dplyr::pull(s_i)
        } else {
          NA_real_
        }
    )

  # join to data_template
  dplyr::left_join(
    data_template,
    res%>%
      dplyr::select(dplyr::any_of(c("id_sample", "s_k_2_l_2", "s_k_2_alpha_2", "s_l_2_alpha_2"))),
    by = "id_sample"
  )

}
