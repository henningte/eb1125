#' Helper function which adds indices to `hpmd_data_stan`
#'
#' Allows to use `tidybayes::recover_type()` with `hpmd_data_stan`, yet
#' avoids recomputing the Stan models when I add a new index (because indices
#' are added only afterwards).
#'
#' @param x The main data table.
#'
#' @param x_stan_data `hpmd_data_stan` (A list).
#'
#' @param ...
#'
#' @export
hpmd_stan_data_add_indices <- function(x, x_stan_data, ...) {

  .dots <- list(...)

  res <- hpmd_prepare_data_stan_1_1(x = x, ...)

  # for bookkeeping
  res <-
    c(
      # same parameters as for the leaching project
      leaching_stan_data_add_indices(x, x_stan_data, ...),
      tibble::lst(
        hpm_microhabitat =
          if(.dots$has_hpm_parameters) {
            res %>%
              dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
              dplyr::pull(hpm_microhabitat2) %>%
              sort()
          } else {
            res %>%
              dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
              dplyr::pull(hpm_microhabitat) %>%
              sort()
          },
        index_samples_to_hpm =
          res %>%
          dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
          dplyr::pull(id_sample_incubation_start)
      )
    )

  res

}


#' Helper function to extract MCMC draws from a single model object
#'
#' @export
hpmd_stan_extract_draws_one_model <- function(x_stan_model, x, x_stan_data, data_template = NULL, id_fit, ...) {

  .dots <- list(...)

  # extract the same parameters as for the leaching project
  res_leaching <- leaching_stan_extract_draws(x_stan_model, x, x_stan_data, data_template = data_template, id_fit, ...)

  if(is.null(data_template)) {
    data_template <- x
  }

  x_stan_data <- hpmd_stan_data_add_indices(x = x, x_stan_data = x_stan_data, ...)

  x_stan_model <- tidybayes::recover_types(x_stan_model, x_stan_data)

  rvar_inv_logit <- posterior::rfun(binomial()$linkinv)

  # parameters on global and hpm sample level
  res_id_sample <-
    list(
      base =
        x_stan_model %>%
        tidybayes::spread_rvars(
          layer_water_table_depth_to_surface_1[index_samples_to_hpm],
          layer_degree_of_saturation_1[index_samples_to_hpm]
        ) %>%
        dplyr::mutate(
          layer_water_table_depth_to_surface_1 = layer_water_table_depth_to_surface_1 * x_stan_data$layer_water_table_depth_to_surface_1_p2 + x_stan_data$layer_water_table_depth_to_surface_1_p1,
          index_samples_to_hpm = x_stan_data$index_samples_to_hpm
        ),
      base_index_test =
        x_stan_model %>%
        tidybayes::spread_rvars(
          hpm_k_2_rep[index_hpm_to_test]
        ) %>%
        dplyr::mutate(
          index_samples_to_hpm = x_stan_data$index_samples_to_hpm[x_stan_data$index_hpm_to_test]
        ) %>%
        dplyr::select(-index_hpm_to_test)
    )

  if(.dots$has_hpm_parameters) {
    res_id_sample$hpm_parameters <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        m68_p1,
        m68_p2,
        m69_p1,
        m69_p2,
        m68_p3_2[index_samples_to_hpm]
      ) %>%
      dplyr::mutate(
        m68_p1 = m68_p1 * x_stan_data$m68_p1_p3,
        m68_p2 = m68_p2 * x_stan_data$m68_p2_p3,
        m69_p2 = m69_p2 * x_stan_data$m69_p2_p3,
        index_samples_to_hpm = x_stan_data$index_samples_to_hpm
      )
  }

  if(.dots$has_hpm_leaching) {
    res_id_sample$hpm_leaching <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        hpm_l_2[index_samples_to_hpm]
      ) %>%
      dplyr::mutate(
        index_samples_to_hpm = x_stan_data$index_samples_to_hpm
      )

    res_id_sample$hpm_leaching_test <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        hpm_l_2_rep[index_hpm_to_test]
      ) %>%
      dplyr::mutate(
        index_samples_to_hpm = x_stan_data$index_samples_to_hpm[x_stan_data$index_hpm_to_test]
      ) %>%
      dplyr::select(-index_hpm_to_test)
  }

  res_id_sample <-
    purrr::reduce(
      res_id_sample,
      dplyr::left_join,
      by = "index_samples_to_hpm"
    ) %>%
    dplyr::mutate(
      index_samples_to_hpm = as.integer(as.character(index_samples_to_hpm))
    ) %>%
    dplyr::rename(
      id_sample_incubation_start = index_samples_to_hpm
    )


  # hpm microhabitats
  res_hpm_microhabitat <- tibble::tibble(id_sample_incubation_start = NA_integer_)
  if(.dots$has_hpm_leaching) {
    res_hpm_microhabitat <-
      x_stan_model %>%
      tidybayes::spread_rvars(
        hpm_l_2_p1[hpm_microhabitat]
      ) %>%
      dplyr::mutate(
        hpm_l_2_p1 = rvar_inv_logit(hpm_l_2_p1 * x_stan_data$hpm_l_2_p1_p2  + x_stan_data$hpm_l_2_p1_p1), #---note: maximum possible initial leaching loss per species at a degree of saturation of 0
      ) %>%
      dplyr::slice(as.integer(x_stan_data$index_hpm_microhabitats_to_hpm)) %>%
      dplyr::mutate(
        id_sample_incubation_start = as.integer(as.character(x_stan_data$index_samples_to_hpm))
      ) %>%
      dplyr::select(-hpm_microhabitat)
  }


  # combine
  res <-
    res_leaching %>%
    dplyr::left_join(
      res_id_sample,
      by = "id_sample_incubation_start"
    ) %>%
    dplyr::left_join(
      res_hpm_microhabitat,
      by = "id_sample_incubation_start"
    )

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
hpmd_stan_extract_draws <- function(x_stan_model, x, x_stan_data, data_template = NULL, id_fit, ...) {

  if(is.data.frame(x_stan_model)) {
    if(nrow(as.data.frame(x_stan_model)) == 0) {#---note: empty model
      return(NULL)
    }
  }

  .dots <- list(...)

  res <-
    if(.dots$has_cross_validation) {
      purrr::map_dfr(seq_along(x_stan_model), function(i) {
        print(i)
        hpmd_stan_extract_draws_one_model(
          x_stan_model = x_stan_model[[i]],
          x = x,
          x_stan_data = x_stan_data[[i]],
          data_template = data_template,
          id_fit = id_fit,
          ...
        ) %>%
          dplyr::mutate(
            id_cross_validation_block_extraction = i
          )
      })
    } else {
      hpmd_stan_extract_draws_one_model(
        x_stan_model = x_stan_model,
        x = x,
        x_stan_data = x_stan_data,
        data_template = data_template,
        id_fit = id_fit,
        ...
      )
    }

  res

}
