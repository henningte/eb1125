# Functions to make predictions with the standard HPM decomposition module and compare them to predictions of model 1-4 from the leaching project

#' Generates predictions with the standard HPM decomposition module
#'
#' @param x A data frame. Corresponds to `leaching_stan_draws_4`.
#'
#' @export
hpmd_get_yhat_fit_1 <- function(x, hpmd_data_hpm_microhabitat, id_fit) {

  # get parameter values
  hpmd_peat_properties <-
    hpmd_get_peat_properties()

  hpmd_hpm_standard_parameter_values <-
    hpmd_get_hpm_standard_parameter_values() %>%
    dplyr::select(-variable_pretty) %>%
    tidyr::pivot_wider(
      values_from = "value",
      names_from = "variable"
    )

  # add information on HPM microhabitat
  res <-
    x %>%
    hpmd_identify_samples_for_hpm() %>%
    hpmd_data_add_hpm_microhabitat(
      hpmd_data_hpm_microhabitat = hpmd_data_hpm_microhabitat
    ) %>%
    hpmd_add_cross_validation_indices()

  # compute peat properties
  res <-
    res %>%
    dplyr::left_join(
      res %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
        dplyr::select(id_sample_incubation_start, sample_depth_upper, water_table_depth_to_surface_average) %>%
        dplyr::mutate(
          layer_total_porosity_1 = {
            mu <- hpmd_peat_properties$layer_total_porosity_1_average
            sigma <- hpmd_peat_properties$layer_total_porosity_1_error
            phi <- ((mu * (1 - mu))/sigma^2 - 1)
            posterior::rvar_rng(rbeta, n = length(id_sample_incubation_start), shape1 = mu * phi, shape2 = (1 - mu) * phi, ndraws = with(leaching_mcmc_config, (niter - nwarmup) * nchains))
          },
          minimum_water_content_at_surface_1 = {
            mu <- hpmd_peat_properties$minimum_water_content_at_surface_1_average
            sigma <- hpmd_peat_properties$minimum_water_content_at_surface_1_error
            phi <- ((mu * (1 - mu))/sigma^2 - 1)
            posterior::rvar_rng(rbeta, n = length(id_sample_incubation_start), shape1 = mu * phi, shape2 = (1 - mu) * phi, ndraws = with(leaching_mcmc_config, (niter - nwarmup) * nchains))
          },
          water_table_depth_to_surface_1 = {
            mu <- water_table_depth_to_surface_average
            sigma <- hpmd_peat_properties$water_table_depth_to_surface_1_error
            posterior::rvar_rng(rnorm, n = length(id_sample_incubation_start), mean = mu, sd = sigma, ndraws = with(leaching_mcmc_config, (niter - nwarmup) * nchains))
            },
          layer_degree_of_saturation_1 =
            m7_rvars(
              layer_depth_midpoint_1 = sample_depth_upper,
              layer_total_porosity_1 = layer_total_porosity_1,
              water_table_depth_to_surface_1 = water_table_depth_to_surface_1,
              minimum_water_content_at_surface_1 = minimum_water_content_at_surface_1,
              layer_total_porosity_at_surface = layer_total_porosity_1
            ) / layer_total_porosity_1,
          layer_water_table_depth_1 = water_table_depth_to_surface_1 - sample_depth_upper,
        ) %>%
        dplyr::select(-sample_depth_upper, -water_table_depth_to_surface_average),
      by = "id_sample_incubation_start"
    )

  # make HPM predictions
  res <-
    res %>%
    dplyr::bind_cols(
      hpmd_hpm_standard_parameter_values %>%
        dplyr::filter(! is.na(m69_p1)) %>%
        dplyr::select(dplyr::all_of(c("m69_p1", "m69_p2", "m68_p1", "m68_p2")))
    ) %>%
    dplyr::left_join(
      hpmd_hpm_standard_parameter_values %>%
        dplyr::filter(! is.na(hpm_microhabitat)) %>%
        dplyr::select(dplyr::all_of(c("hpm_microhabitat", "m68_p3"))),
      by = "hpm_microhabitat"
    ) %>%
    dplyr::mutate(
      hpm_k_2_rep =
        m68(
          layer_degree_of_saturation_1 = layer_degree_of_saturation_1,
          layer_water_table_depth_1 = layer_water_table_depth_1,
          m68_p3 = m68_p3,
          m68_p4 = 1,
          m68_p5 = 1
        )
    )

  # add new id_fit
  res <-
    res %>%
    dplyr::mutate(
      id_fit = !!id_fit
    )

  # export
  res_file <- "_targets_rvars/hpmd_yhat_1.rds"
  saveRDS_rvars(res, res_file)

  res_file

}
