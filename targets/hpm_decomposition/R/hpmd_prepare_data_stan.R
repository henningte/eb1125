#' Helper function to define factors where necessary
#'
#' @param x The main data table
#'
#' @param ...
hpmd_prepare_data_stan_1_1 <- function(x, ...) {

  .dots <- list(...)

  if(! .dots$has_outlier) {
    x <-
      dplyr::left_join(
        x,
        leaching_identify_outliers_mass_relative_mass_mu(
          x_stan_draws = readRDS_rvars(.dots$hpmd_stan_draws_1),
          alpha = 0.99
        ),
        by = "id_sample_incubation_start"
      ) %>%
      dplyr::mutate(
        index_hpm = index_hpm & ! is_outlier
      ) %>%
      dplyr::select(-is_outlier)
  }

  res <- leaching_prepare_data_stan_1_1(x, ...)

  # define indices
  res <-
    res %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("hpm_microhabitat", "hpm_microhabitat2")),
        function(.x) {
          factor(.x, levels = unique(.x[index_hpm]))
        }
      )
    )

  res

}

#' Prepares data for Stan models
#'
#' @param x The main data table
#'
#' @param ...
#'
#' @param mass_relative_mass_offset
#'
#' @export
hpmd_prepare_data_stan <- function(x, ..., mass_relative_mass_offset) {

  .dots <- list(...)

  # get parameter values
  hpmd_peat_properties <-
    hpmd_get_peat_properties() %>%
    dplyr::mutate(
      layer_total_porosity_1_phi = {
        mu <- layer_total_porosity_1_average
        sigma <- layer_total_porosity_1_error
        ((mu * (1 - mu))/sigma^2 - 1)
      },
      minimum_water_content_at_surface_1_phi = {
        mu <- minimum_water_content_at_surface_1_average
        sigma <- minimum_water_content_at_surface_1_error
        ((mu * (1 - mu))/sigma^2 - 1)
      }
    )

  hpmd_hpm_standard_parameter_values <-
    hpmd_get_hpm_standard_parameter_values() %>%
    dplyr::select(-variable_pretty) %>%
    tidyr::pivot_wider(
      values_from = "value",
      names_from = "variable"
    )

  # add factors
  .x <- hpmd_prepare_data_stan_1_1(x = x, ...)

  # add constants
  .x <-
    .x  %>%
    dplyr::left_join(
      hpmd_hpm_standard_parameter_values %>%
        dplyr::filter(! is.na(hpm_microhabitat)) %>%
        dplyr::select(dplyr::all_of(c("hpm_microhabitat", "m68_p3"))),
      by = "hpm_microhabitat"
    )

  # parameters from this project
  res <-
    tibble::lst(
      # controls
      has_hpm_parameters = .dots$has_hpm_parameters,
      has_hpm_leaching = .dots$has_hpm_leaching,
      has_cross_validation = .dots$has_cross_validation,
      # counter
      N_hpm =
        .x %>%
        dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
        nrow(),
      N_hpm_microhabitats =
        if(.dots$has_hpm_parameters) {
          .x %>%
            dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
            dplyr::pull(hpm_microhabitat2) %>%
            unique() %>%
            length()
        } else {
          .x %>%
            dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
            dplyr::pull(hpm_microhabitat) %>%
            unique() %>%
            length()
        },
      N_train = N_hpm,
      N_test = N_hpm,
      # indices
      index_samples_to_hpm =
        .x %>%
        dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
        dplyr::pull(id_sample_incubation_start),
      index_hpm_microhabitats_to_hpm =
        if(.dots$has_hpm_parameters) {
          .x %>%
            dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
            dplyr::pull(hpm_microhabitat2)
        } else {
          .x %>%
            dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
            dplyr::pull(hpm_microhabitat)
        },
      index_hpm_to_train =
        .x %>%
        dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
        nrow() %>%
        seq_len(),
      index_hpm_to_test = index_hpm_to_train,
      # peat properties
      layer_total_porosity_1_p1 = hpmd_peat_properties$layer_total_porosity_1_average * hpmd_peat_properties$layer_total_porosity_1_phi, # beta
      layer_total_porosity_1_p2 = (1 - hpmd_peat_properties$layer_total_porosity_1_average) * hpmd_peat_properties$layer_total_porosity_1_phi,
      layer_water_table_depth_to_surface_1_p1 =
        .x %>%
        dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
        dplyr::pull(water_table_depth_to_surface_average), # normal
      layer_water_table_depth_to_surface_1_p2 =
        hpmd_peat_properties$water_table_depth_to_surface_1_error,
      layer_minimum_degree_of_saturation_at_surface_1_p1 = hpmd_peat_properties$minimum_water_content_at_surface_1_average * hpmd_peat_properties$minimum_water_content_at_surface_1_phi, # beta
      layer_minimum_degree_of_saturation_at_surface_1_p2 = (1 - hpmd_peat_properties$minimum_water_content_at_surface_1_average) * hpmd_peat_properties$minimum_water_content_at_surface_1_phi,
      layer_depth_midpoint_1 =
        .x %>%
        dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
        dplyr::pull(sample_depth_upper),
      hpm_k_2_p1_p1 = 20, # 10, # gamma
      hpm_k_2_p1_p2 = 20/20, # 10/10,
      hpm_k_2_p1_p3 = 20, # 1,
      hpm_k_2_p3 = 20,
      # hpm parameters (dummy values)
      m69_p1_p1 = 1, # beta
      m69_p1_p2 = 1,
      m69_p2_p1 = 1, # gamma
      m69_p2_p2 = 1,
      m69_p2_p3 = 1,
      m68_p1_p1 = 1, # gamma
      m68_p1_p2 = 1,
      m68_p1_p3 = 1,
      m68_p2_p1 = 1, # gamma
      m68_p2_p2 = 1,
      m68_p2_p3 = 1,
      m68_p3_2_p1_p1 = 1, # normal
      m68_p3_2_p1_p2 = 1,
      m68_p3_2_p2_p1 = 1, # 0
      # group-level standard deviations (dummy values)
      m68_p3_2_p1_p2_p2 = 1,
      m68_p3_2_p2_p2_p2 = 1, # 0.2
      # hpm leaching (dummy values)
      hpm_l_2_p4_p1 = 1,
      hpm_l_2_p4_p2 = 1,
      hpm_l_2_p4_p3 = 1,
      hpm_l_2_p1_p1 = 1,
      hpm_l_2_p1_p2 = 1,
      hpm_l_2_p2_p1 = 1,
      hpm_l_2_p3_p1 = 1,
      hpm_l_2_p3_p2 = 1,
      # group-level standard deviations (dummy values)
      hpm_l_2_p1_p2_p2 = 1,
      hpm_l_2_p2_p2_p2 = 1,
      # constants
      m69_p1_constant = hpmd_hpm_standard_parameter_values$m69_p1[[1]],
      m69_p2_constant = hpmd_hpm_standard_parameter_values$m69_p2[[1]],
      m68_p1_constant = hpmd_hpm_standard_parameter_values$m68_p1[[1]],
      m68_p2_constant = hpmd_hpm_standard_parameter_values$m68_p2[[1]],
      m68_p3_2_constant =
        .x %>%
        dplyr::filter(! duplicated(id_sample_incubation_start) & index_hpm) %>%
        dplyr::pull(m68_p3)
    )

  # add HPM parameters
  if(.dots$has_hpm_parameters) {
    # hpm parameters
    res$m69_p1_p1 = hpmd_hpm_standard_parameter_values$m69_p1[[1]] * 30 # beta
    res$m69_p1_p2 = (1 - hpmd_hpm_standard_parameter_values$m69_p1[[1]]) * 30
    res$m69_p2_p1 = 20 # gamma
    res$m69_p2_p2 = 20/hpmd_hpm_standard_parameter_values$m69_p2[[1]]
    res$m69_p2_p3 = 5
    res$m68_p1_p1 = 5 # gamma
    res$m68_p1_p2 = 5/hpmd_hpm_standard_parameter_values$m68_p1[[1]]
    res$m68_p1_p3 = 1/200
    res$m68_p2_p1 = 5 # gamma
    res$m68_p2_p2 = 5/hpmd_hpm_standard_parameter_values$m68_p2[[1]]
    res$m68_p2_p3 = 1
    res$m68_p3_2_p1_p1 = -2.2 # normal
    res$m68_p3_2_p1_p2 = 0.3
    res$m68_p3_2_p2_p1 = 0
    # group-level standard deviations
    res$m68_p3_2_p1_p2_p2 = 0.3
    res$m68_p3_2_p2_p2_p2 = 0.3 # 0.2
  }

  # add HPM leaching parameters
  if(.dots$has_hpm_leaching) {
    # hpm leaching
    res$hpm_l_2_p4_p1 = 10
    res$hpm_l_2_p4_p2 = 10/40
    res$hpm_l_2_p4_p3 = 10
    res$hpm_l_2_p1_p1 = -2.2
    res$hpm_l_2_p1_p2 = 0.3
    res$hpm_l_2_p2_p1 = 0
    res$hpm_l_2_p3_p1 = 0.0
    res$hpm_l_2_p3_p2 = 0.5
    # group-level standard deviations
    res$hpm_l_2_p1_p2_p2 = 0.5
    res$hpm_l_2_p2_p2_p2 = 0.5
  }

  # add parameters from this project
  res <-
    c(
      leaching_prepare_data_stan(
        x = .x,
        ...,
        mass_relative_mass_offset = mass_relative_mass_offset
      ), # same data as for the leaching project
      res
    )

  # for cross-validation
  if(.dots$has_cross_validation) {

    res <-
      purrr::map(unique(.x$hpm_cross_validation_id_block[.x$index_hpm_cross_validation]), function(i) {

        res$N_train <-
          .x %>%
          dplyr::filter(index_hpm & (hpm_cross_validation_id_block != i | is.na(hpm_cross_validation_id_block))) %>%
          dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
          nrow()
        res$N_test <-
          .x %>%
          dplyr::filter(index_hpm & hpm_cross_validation_id_block == i) %>%
          dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
          nrow()
        res$index_hpm_to_train <-
          .x %>%
          dplyr::left_join(
            .x %>%
              dplyr::filter(index_hpm) %>%
              dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
              dplyr::mutate(
                index = seq_len(nrow(.))
              ) %>%
              dplyr::select(id_sample, index),
            by = "id_sample"
          ) %>%
          dplyr::filter(index_hpm & (hpm_cross_validation_id_block != i | is.na(hpm_cross_validation_id_block))) %>%
          dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
          dplyr::pull(index)
        res$index_hpm_to_test <-
          .x %>%
          dplyr::left_join(
            .x %>%
              dplyr::filter(index_hpm) %>%
              dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
              dplyr::mutate(
                index = seq_len(nrow(.))
              ) %>%
              dplyr::select(id_sample, index),
            by = "id_sample"
          ) %>%
          dplyr::filter(index_hpm & hpm_cross_validation_id_block == i) %>%
          dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
          dplyr::pull(index)

        if(length(res$index_hpm_to_test) == 1) {
          res$index_hpm_to_test <- array(res$index_hpm_to_test, dim = 1)
        }

        res

      })

  }

  res

}
