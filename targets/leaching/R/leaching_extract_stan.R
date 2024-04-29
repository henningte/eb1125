#' Helper function which adds indices to `leaching_data_stan`
#'
#' Allows to use `tidybayes::recover_type()` with `leaching_data_stan`, yet
#' avoids recomputing the Stan models when I add a new index (because indices
#' are added only afterwards).
#'
#' @param x The main data table.
#'
#' @param x_stan_data `leaching_data_stan` (A list).
#'
#' @param ...
#'
#' @export
leaching_stan_data_add_indices <- function(x, x_stan_data, ...) {

  .dots <- list(...)

  res <- leaching_prepare_data_stan_1_1(x = x, ...)

  # for bookkeeping
  res <-
    c(
      x_stan_data,
      tibble::lst(
        id_sample =
          res %>%
          dplyr::pull(id_sample),
        id_sample_incubation_start =
          res %>%
          dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
          dplyr::pull(id_sample_incubation_start),
        taxon_rank_value =
          res %>%
          dplyr::filter(! duplicated(taxon_rank_value)) %>%
          dplyr::pull(taxon_rank_value) %>%
          sort(),
        speciesxstudies =
          res %>%
          dplyr::filter(! duplicated(speciesxstudies)) %>%
          dplyr::pull(speciesxstudies) %>%
          sort()
      )
    )

  res

}




#' Extracts parameters from fitted Stan models
#'
#' @param x_stan_model A stanfit object.
#'
#' @param x_stan_data The Stan data used to compute `x_stan_model`.
#'
#' @param x The main data table.
#'
#' @param data_template A data frame with columns needed to bind the extracted
#' parameters. (Here, this is mainly `id_sample`, `id_sample_incubation_start`).
#'
#' @param id_fit Integer. The current ID of `x_stan_model`.
#'
#' @param ... Further arguments, options from `leaching_d_models`.
#'
#' @export
leaching_stan_extract_draws <- function(x_stan_model, x, x_stan_data, data_template = NULL, id_fit, ...) {

  .dots <- list(...)

  if(is.null(data_template)) {
    data_template <- x
  }

  x_stan_data <- leaching_stan_data_add_indices(x = x, x_stan_data = x_stan_data, ...)

  x_stan_model <- tidybayes::recover_types(x_stan_model, x_stan_data)

  rvar_inv_logit <- posterior::rfun(binomial()$linkinv)

  # parameters on sample level
  res_id_sample <-
    x_stan_model %>%
    tidybayes::spread_rvars(
      m_rep[id_sample],
      mu[id_sample],
      phi[id_sample]
    ) %>%
    dplyr::mutate(
      id_sample = x_stan_data$id_sample[id_sample]
    )

  # parameters on experiment level
  res_id_sample_incubation_start <-
    list(
      base =
        x_stan_model %>%
        tidybayes::spread_rvars(
          k_2[id_sample_incubation_start],
          phi_2_p2[id_sample_incubation_start]
        )
    )

  if(.dots$has_leaching) {
    res_id_sample_incubation_start$l_2 <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        l_2[id_sample_incubation_start]
      )
  }

  if(.dots$has_alpha) {
    res_id_sample_incubation_start$alpha_2 <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        alpha_2[id_sample_incubation_start]
      )
  }

  res_id_sample_incubation_start <-
    purrr::reduce(
      res_id_sample_incubation_start,
      dplyr::left_join,
      by = "id_sample_incubation_start"
    ) %>%
    dplyr::mutate(
      id_sample_incubation_start = as.integer(as.character(x_stan_data$id_sample_incubation_start[id_sample_incubation_start]))
    )

  # parameters on species level
  res_taxon_rank_value <-
    list(
      base =
        x_stan_model %>%
        tidybayes::spread_rvars(
          k_2_p1,
          k_2_p1_p2,
          k_2_p2[taxon_rank_value],
          k_2_p2_p2[taxon_rank_value]
        ) %>%
        dplyr::mutate(
          k_2_p1_p2 = k_2_p1_p2 * x_stan_data$k_2_p1_p2_p1,
          k_2_p2_p2 = k_2_p2_p2 * x_stan_data$k_2_p2_p2_p1,
          k_2_species = exp(k_2_p1 * k_2_p1_p2 + x_stan_data$k_2_p1_p1 + k_2_p2 * k_2_p2_p2 + x_stan_data$k_2_p2_p1)

        )
    )

  if(.dots$has_leaching) {
    res_taxon_rank_value$l_2 <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        l_2_p1,
        l_2_p1_p2,
        l_2_p2[taxon_rank_value],
        l_2_p2_p2[taxon_rank_value]
      ) %>%
      dplyr::mutate(
        l_2_p1_p2 = l_2_p1_p2 * x_stan_data$l_2_p1_p2_p1,
        l_2_p2_p2 = l_2_p2_p2 * x_stan_data$l_2_p2_p2_p1,
        l_2_species = rvar_inv_logit(l_2_p1 * l_2_p1_p2 + x_stan_data$l_2_p1_p1 + l_2_p2 * l_2_p2_p2 + x_stan_data$l_2_p2_p1)
      )
  }

  if(.dots$has_alpha) {
    res_taxon_rank_value$alpha_2 <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        alpha_2_p1,
        alpha_2_p2[taxon_rank_value]
      ) %>%
      dplyr::mutate(
        alpha_2_species = 1 + exp(alpha_2_p1 * x_stan_data$alpha_2_p1_p2 + x_stan_data$alpha_2_p1_p1 + alpha_2_p2 * x_stan_data$alpha_2_p2_p2 + x_stan_data$alpha_2_p2_p1)
      ) %>%
      dplyr::select(-alpha_2_p1, -alpha_2_p2)
  }

  res_taxon_rank_value <-
    purrr::reduce(
      res_taxon_rank_value,
      dplyr::left_join,
      by = "taxon_rank_value"
    ) %>%
    dplyr::mutate(
      taxon_rank_value = as.character(x_stan_data$taxon_rank_value)[taxon_rank_value]
    )


  # parameters on speciesxstudies level
  res_speciesxstudies <-
    list(
      base =
        x_stan_model %>%
        tidybayes::spread_rvars(
          k_2_p3_p2[speciesxstudies]
        ) %>%
        dplyr::mutate(
          k_2_p3_p2 = k_2_p3_p2 * x_stan_data$k_2_p3_p2_p1
        )
    )

  if(.dots$has_leaching) {
    res_speciesxstudies$l_2 <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        l_2_p3_p2[speciesxstudies]
      ) %>%
      dplyr::mutate(
        l_2_p3_p2 = l_2_p3_p2 * x_stan_data$l_2_p3_p2_p1
      )
  }

  res_speciesxstudies <-
    purrr::reduce(
      res_speciesxstudies,
      dplyr::left_join,
      by = "speciesxstudies"
    ) %>%
    dplyr::mutate(
      speciesxstudies = as.character(x_stan_data$speciesxstudies)[speciesxstudies]
    )

  # combine
  res <-
    data_template %>%
    dplyr::left_join(
      res_id_sample,
      by = "id_sample"
    ) %>%
    dplyr::left_join(
      res_id_sample_incubation_start,
      by = "id_sample_incubation_start"
    ) %>%
    dplyr::left_join(
      res_taxon_rank_value,
      by = "taxon_rank_value"
    ) %>%
    dplyr::left_join(
      res_speciesxstudies,
      by = "speciesxstudies"
    )

  # add current id_fit
  res <-
    res %>%
    dplyr::mutate(
      id_fit = id_fit
    )

  res

}



#' Combines extracted values with another data frame and extracts columns
#'
#' Both data frames are joined via `c("id_sample", "id_sample_incubation_start", "taxon_rank_value", "speciesxstudies")`.
#'
#' @param x
#'
#' @param x_stan_draws A data frame with Stan draws.
#'
#' @param ... Arguments passed to `dplyr::select()` after the join.
#'
#' @export
leaching_combine_x_with_stan_draws <- function(x, x_stan_draws, ...) {

  .dots <- list(...)

  # join
  res <-
    dplyr::left_join(
      x,
      x_stan_draws,
      by = c("id_sample", "id_sample_incubation_start", "taxon_rank_value", "speciesxstudies")
    )

  if(length(.dots) > 0) {
    res <-
      res %>%
      dplyr::select(...)
  }

  res

}


