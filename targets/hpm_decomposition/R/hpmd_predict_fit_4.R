#' Predicts decomposition rates and initial leaching losses for fit 4 (model 3-4)
#'
#' @examples
#' newdata <-
#'   tibble::tibble(
#'     incubation_duration = seq(0, 100, 1),
#'     m0 = 1,
#'     layer_degree_of_saturation_1 = 0.4,
#'     layer_water_table_depth_to_surface_1 = 20,
#'     sample_depth_lower = 10,
#'     hpm_taxon_rank_value = "Sphagnum fuscum"
#'  )
#'
#' @export
hpmd_predict_fit_4 <- function(newdata, x_stan_fit_4, x_stan_data_4, hpmd_data) {

  m <- x_stan_fit_4
  m_config <- x_stan_data_4

  rvar_inv_logit <- posterior::rfun(binomial()$linkinv)

  # transformed data
  newdata <-
    newdata %>%
    dplyr::mutate(
      layer_water_table_depth_1 = layer_water_table_depth_to_surface_1 - sample_depth_lower
    )

  # extract parameters
  m_parameters_alpha <-
    m %>%
    tidybayes::spread_rvars(
      alpha_2_p1,
      alpha_2_p2[index_taxon_rank_value]
    ) %>%
    dplyr::mutate(
      alpha_2_p1 = alpha_2_p1 * m_config$alpha_2_p1_p2 + m_config$alpha_2_p1_p1,
      alpha_2_p2 = alpha_2_p2 * m_config$alpha_2_p2_p2 + m_config$alpha_2_p2_p1,
      alpha_2 = 1.0 + exp(alpha_2_p1 + alpha_2_p2) #---todo: how to predict for new studies and samples here?
    ) %>%
    dplyr::slice(as.integer(m_config$index_species_to_samples)) %>%
    dplyr::slice(as.integer(m_config$index_samples_to_measurements)) %>%
    dplyr::mutate(
      hpm_microhabitat2 = hpmd_data$hpm_microhabitat2[hpmd_data$id_citation != "Bengtsson.2017"]
    ) %>%
    dplyr::filter(! duplicated(hpm_microhabitat2) & ! is.na(hpm_microhabitat2))

  m_parameters_phi <-
    m %>%
    tidybayes::spread_rvars(
      phi_2_p2_p1,
      phi_2_p2_p2[index_taxon_rank_value],
      phi_2_p1
    ) %>%
    dplyr::mutate(
      phi_2_p2_p1 = phi_2_p2_p1 * m_config$phi_2_p2_p1_p2 + m_config$phi_2_p2_p1_p1,
      phi_2_p2_p2 = phi_2_p2_p2 * m_config$phi_2_p2_p2_p2 + m_config$phi_2_p2_p2_p1,
      phi_2_p2 = exp(phi_2_p2_p1 + phi_2_p2_p2), #---todo: how to predict for new studies and samples here?
      phi_2_p1 = phi_2_p1 * m_config$phi_2_p1_p3,
      phi_2 = posterior::rvar_rng(rgamma, n = length(phi_2_p2), shape = phi_2_p1, rate = phi_2_p1/phi_2_p2)
    ) %>%
    dplyr::slice(as.integer(m_config$index_species_to_samples)) %>%
    dplyr::slice(as.integer(m_config$index_samples_to_measurements)) %>%
    dplyr::mutate(
      hpm_microhabitat2 = hpmd_data$hpm_microhabitat2[hpmd_data$id_citation != "Bengtsson.2017"]
    ) %>%
    dplyr::filter(! duplicated(hpm_microhabitat2) & ! is.na(hpm_microhabitat2))

  m_parameters <-
    m %>%
    tidybayes::spread_rvars(
      m69_p1,
      m69_p2,
      m68_p1,
      m68_p2,
      m68_p3_2_p1[index_hpm_taxon_rank_value],
      hpm_l_2_p3,
      hpm_l_2_p1[index_hpm_taxon_rank_value],
      hpm_k_2_p1,
      hpm_l_2_p4
    ) %>%
    dplyr::mutate(
      m69_p2 = m69_p2 * m_config$m69_p2_p3,
      m68_p1 = m68_p1 * m_config$m68_p1_p3,
      m68_p2 = m68_p2 * m_config$m68_p2_p3,
      m68_p3_2 = exp(m68_p3_2_p1 * m_config$m68_p3_2_p1_p2 + m_config$m68_p3_2_p1_p1),
      hpm_l_2_p3 = hpm_l_2_p3 * m_config$hpm_l_2_p3_p2 + m_config$hpm_l_2_p3_p1,
      hpm_l_2_p1 = hpm_l_2_p1 * m_config$hpm_l_2_p1_p2 + m_config$hpm_l_2_p1_p1,
      hpm_l_2_p4 = hpm_l_2_p4 * m_config$hpm_l_2_p4_p3,
      hpm_k_2_p1 = hpm_k_2_p1 * m_config$hpm_k_2_p1_p3,
      hpm_microhabitat2 = levels(m_config$index_hpm_microhabitats_to_hpm)
    ) %>%
    dplyr::left_join(
      m_parameters_alpha %>%
        dplyr::select(alpha_2, hpm_microhabitat2),
      by = "hpm_microhabitat2"
    ) %>%
    dplyr::left_join(
      m_parameters_phi %>%
        dplyr::select(phi_2, hpm_microhabitat2),
      by = "hpm_microhabitat2"
    )

  # predict decomposition rate
  newdata <-
    newdata %>%
    dplyr::mutate(
      hpm_microhabitat2 =
        factor(hpm_taxon_rank_value, levels = levels(m_config$index_hpm_microhabitats_to_hpm))
    ) %>%
    dplyr::left_join(
      m_parameters %>%
        dplyr::select(! dplyr::all_of(setdiff(intersect(colnames(newdata), colnames(m_parameters)), "hpm_microhabitat2"))),
      by = "hpm_microhabitat2"
    ) %>%
    dplyr::mutate(
      hpm_k_2 =
        m68(
          layer_degree_of_saturation_1 = newdata$layer_degree_of_saturation_1,
          layer_water_table_depth_1 = newdata$layer_water_table_depth_1,
          m69_p1 = m69_p1,
          m69_p2 = m69_p2,
          m68_p1 = m68_p1,
          m68_p2 = m68_p2,
          m68_p3 = m68_p3_2,
          m68_p4 = 1,
          m68_p5 = 1
        ),
      k0 =
        posterior::rvar_rng(
          rgamma,
          n = length(hpm_k_2_p1),
          shape = hpm_k_2_p1,
          rate = hpm_k_2_p1 / hpm_k_2
        )
    )

  # predict initial leaching loss
  newdata <-
    newdata %>%
    dplyr::mutate(
      l0 =
        if("l0" %in% colnames(.)) {
          posterior::as_rvar(l0)
        } else {
          mu <- rvar_inv_logit(hpm_l_2_p1 + hpm_l_2_p3 * layer_degree_of_saturation_1)
          posterior::rvar_rng(
            rbeta,
            n = length(mu),
            shape1 = mu * hpm_l_2_p4,
            shape2 = (1 - mu) * hpm_l_2_p4
          )
        }
    )

  # predict remaining mass
  newdata <-
    newdata %>%
    dplyr::mutate(
      l0 =
        dplyr::case_when(
          incubation_duration <= 0 ~ posterior::as_rvar(0),
          TRUE ~ .data$l0
        ),
      mass_relative_mass_mu =
        (m0 - .data$l0) / (1.0 + (.data$alpha_2 - 1.0) * .data$k0 * .data$incubation_duration)^(1.0 / (.data$alpha_2 - 1.0)),
      mass_relative_mass_mu =
        posterior::rvar_ifelse(mass_relative_mass_mu >= 1.0, 1.0 - m_config$s, mass_relative_mass_mu),
      mass_relative_mass =
        posterior::rvar_rng(
          rbeta,
          n = length(mass_relative_mass_mu),
          shape1 = mass_relative_mass_mu * phi_2,
          shape2 = (1 - mass_relative_mass_mu) * phi_2
        )
    )

  newdata

}
