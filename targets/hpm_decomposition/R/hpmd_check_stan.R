#### Combine all ####

#' Wrapper function which combines all checks
#'
#' @export
hpmd_stan_get_all_checks <- function(
    prior_models,
    posterior_models,
    x_stan_models
) {

  prior_models <- purrr::map(prior_models, readRDS_rvars)
  posterior_models <- purrr::map(posterior_models, readRDS_rvars)

  leaching_d_models <- hpmd_d_models

  res <-
    tibble::lst(
      plots =
        list(
          # m_rep: prior
          hpmd_plot_ppc_prior_m =
            if(length(prior_models) > 0) {
              prior_models %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "mass_relative_mass",
                  yrep_name = "m_rep",
                  x_lab = "Remaining mass (mass %)",
                  x_scale = 100,
                  file = "figures/hpmd_plot_ppc_prior_m.pdf",
                  width = 10,
                  height = 3,
                  add_prefix_model = FALSE
                )
            } else {
              NULL
            },
          # phi: prior
          hpmd_plot_ppc_prior_phi =
            if(length(prior_models) > 0) {
              prior_models %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "mass_relative_mass_precision",
                  yrep_name = "phi",
                  x_lab = "&Phi; (-)",
                  x_scale = 1,
                  file = "figures/hpmd_plot_ppc_prior_phi.pdf",
                  width = 10,
                  height = 3,
                  add_prefix_model = FALSE,
		  x_axis_use_log_scale = TRUE,
		  x_axis_label_digits = 0
                )
            } else {
              NULL
            },
          # hpm_k_2_rep: prior
          hpmd_plot_ppc_prior_hpm_k_2 =
            {
              hpmd_stan_draws_1 <- readRDS_rvars(tar_read(hpmd_stan_draws_1))

              list(
                readRDS_rvars(tar_read(hpmd_stan_draws_6)),
                readRDS_rvars(tar_read(hpmd_stan_draws_7)),
                readRDS_rvars(tar_read(hpmd_stan_draws_8))
              ) %>%
                purrr::map(function(.x) {
                  .x %>%
                    dplyr::mutate(
                      k_2_leaching_id_fit_4 = mean(hpmd_stan_draws_1$k_2)
                    ) %>%
                    dplyr::filter(index_hpm)
                }) %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "k_2_leaching_id_fit_4",
                  yrep_name = "hpm_k_2_rep",
                  x_lab = "<i>k<sub>0</sub></i> (yr<sup>-1</sup>)",
                  x_scale = 1,
                  file = "figures/hpmd_plot_ppc_prior_hpm_k_2.pdf",
                  width = 10,
                  height = 3,
                  add_prefix_model = FALSE,
                  x_axis_label_digits = 2
                )
            },
          # hpm_l_2_rep: prior
          hpmd_plot_ppc_prior_hpm_l_2 =
            {
              hpmd_stan_draws_1 <- readRDS_rvars(tar_read(hpmd_stan_draws_1))
              list(
                readRDS_rvars(tar_read(hpmd_stan_draws_8))
              ) %>%
                purrr::map(function(.x) {
                  .x %>%
                    dplyr::mutate(
                      l_2_leaching_id_fit_4 = mean(hpmd_stan_draws_1$l_2)
                    ) %>%
                    dplyr::filter(index_hpm)
                }) %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "l_2_leaching_id_fit_4",
                  yrep_name = "hpm_l_2_rep",
                  x_lab = "<i>l<sub>0</sub></i> (mass %)",
                  x_scale = 100,
                  file = "figures/hpmd_plot_ppc_prior_hpm_l_2.pdf",
                  width = 10,
                  height = 3,
                  add_prefix_model = FALSE
                )
            },
          # m_rep: posterior
          hpmd_plot_ppc_posterior_m =
            {
              leaching_d_models <- hpmd_d_models
              posterior_models %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "mass_relative_mass",
                  yrep_name = "m_rep",
                  x_lab = "Remaining mass (mass %)",
                  x_scale = 100,
                  file = "figures/hpmd_plot_ppc_posterior_m.pdf",
                  width = 10,
                  height = 6,
                  add_prefix_model = FALSE
                )
            },
          # phi: posterior
          hpmd_plot_ppc_posterior_phi =
            {
              leaching_d_models <- hpmd_d_models
              posterior_models %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "mass_relative_mass_precision",
                  yrep_name = "phi",
                  x_lab = "&Phi; (-)",
                  x_scale = 1,
                  file = "figures/hpmd_plot_ppc_posterior_phi.pdf",
                  width = 10,
                  height = 6,
                  add_prefix_model = FALSE,
		  x_axis_use_log_scale = TRUE,
		  x_axis_label_digits = 0
                )
            },
          # hpm_k_2_rep: posterior
          hpmd_plot_ppc_posterior_hpm_k_2 =
            {
              posterior_models[-1] %>%
                purrr::map(function(.x) {
                  .x %>%
                    dplyr::mutate(
                      k_2 = mean(k_2)
                    ) %>%
                    dplyr::filter(index_hpm)
                }) %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "k_2",
                  yrep_name = "hpm_k_2_rep",
                  x_lab = "<i>k<sub>0</sub></i> (yr<sup>-1</sup>)",
                  x_scale = 1,
                  file = "figures/hpmd_plot_ppc_posterior_hpm_k_2.pdf",
                  width = 12,
                  height = 3,
                  add_prefix_model = FALSE,
                  x_axis_label_digits = 2
                )
            },
          # hpm_l_2_rep: posterior
          hpmd_plot_ppc_posterior_hpm_l_2 =
            {
              posterior_models[c(4, 5)] %>%
                purrr::map(function(.x) {
                  .x %>%
                    dplyr::mutate(
                      l_2 = mean(l_2)
                    ) %>%
                    dplyr::filter(index_hpm)
                }) %>%
                leaching_check_stan_plot_ppc_combined(
                  y_name = "l_2",
                  yrep_name = "hpm_l_2_rep",
                  x_lab = "<i>l<sub>0</sub></i> (mass-%)",
                  x_scale = 100,
                  file = "figures/hpmd_plot_ppc_posterior_hpm_l_2.pdf",
                  width = 10,
                  height = 3,
                  add_prefix_model = FALSE
                )
            },
          # mrep: y vs yrep
          hpmd_plot_y_vs_yrep_m =
            {
              leaching_d_models <- hpmd_d_models
              posterior_models %>%
                leaching_check_stan_plot_y_vs_yrep_combined(
                  y_name = "mass_relative_mass",
                  yrep_name = "m_rep",
                  x_scale = 100,
                  axis_title_unit = "mass %",
                  use_coord_fixed = FALSE,
                  file = "figures/hpmd_plot_y_vs_yrep_m.pdf",
                  width = 11,
                  height = 7
                )
            },
          # phi: y vs yrep
          hpmd_plot_y_vs_yrep_phi =
            {
              leaching_d_models <- hpmd_d_models
              posterior_models %>%
                leaching_check_stan_plot_y_vs_yrep_combined(
                  y_name = "mass_relative_mass_precision",
                  yrep_name = "phi",
                  x_scale = 1,
                  axis_title_unit = "-",
                  use_coord_fixed = FALSE,
                  file = "figures/hpmd_plot_y_vs_yrep_phi.pdf",
                  width = 11,
                  height = 7
                )
            },
          # hpm_k_2_rep: y vs yrep
          hpmd_plot_y_vs_yrep_hpm_k_2 =
            {
              posterior_models[-1] %>%
                purrr::map(function(.x) {
                  .x %>%
                    dplyr::mutate(
                      k_2 = mean(k_2)
                    ) %>%
                    dplyr::filter(index_hpm)
                }) %>%
                leaching_check_stan_plot_y_vs_yrep_combined(
                  y_name = "k_2",
                  yrep_name = "hpm_k_2_rep",
                  x_scale = 1,
                  axis_title_unit = "yr<sup>-1</sup>",
                  use_coord_fixed = FALSE,
                  file = "figures/hpmd_plot_y_vs_yrep_hpm_k_2.pdf",
                  width = 11,
                  height = 3.5
                )
            },
          # hpm_k_2_rep: y vs yrep
          hpmd_plot_y_vs_yrep_hpm_l_2 =
            {
              posterior_models[c(4, 5)] %>%
                purrr::map(function(.x) {
                  .x %>%
                    dplyr::mutate(
                      l_2 = mean(l_2)
                    ) %>%
                    dplyr::filter(index_hpm)
                }) %>%
                leaching_check_stan_plot_y_vs_yrep_combined(
                  y_name = "l_2",
                  yrep_name = "hpm_l_2_rep",
                  x_scale = 100,
                  axis_title_unit = "mass %",
                  use_coord_fixed = FALSE,
                  file = "figures/hpmd_plot_y_vs_yrep_hpm_l_2.pdf",
                  width = 11,
                  height = 3.5
                )
            }
        ),
      tables =
        tibble::lst(
          hpmd_table_mcmc_diagnostics =
            purrr::map2(x_stan_models, hpmd_d_models$id_fit[! hpmd_d_models$prior_only & ! hpmd_d_models$has_cross_validation], function(.x, .y) {
              leaching_check_stan_diagnostics(
                x_stan_model = .x,
                id_fit = .y
              )
            }),
          hpmd_table_mcmc_diagnostics_summary =
            purrr::map2_dfr(x_stan_models, seq_along(hpmd_d_models$id_fit[! hpmd_d_models$prior_only & ! hpmd_d_models$has_cross_validation]), function(.x, .y) {
              leaching_check_stan_diagnostics_summary(
                x_stan_model = .x,
                leaching_table_mcmc_diagnostics = hpmd_table_mcmc_diagnostics[[.y]],
                id_fit = hpmd_d_models$id_fit[! hpmd_d_models$prior_only & ! hpmd_d_models$has_cross_validation][[.y]]
              )
            }),
          hpmd_table_mcse_summary =
            purrr::map_dfr(hpmd_table_mcmc_diagnostics, function(.x) {

              .x %>%
                dplyr::filter(stringr::str_detect(par, "^l_2\\[") | stringr::str_detect(par, "^k_2\\[") | stringr::str_detect(par, "^alpha_2\\[") | stringr::str_detect(par, "^m_rep\\[") | stringr::str_detect(par, "^hpm_k_2\\[") | stringr::str_detect(par, "^hpm_l_2\\[") | stringr::str_detect(par, "^m68_p1\\[") | stringr::str_detect(par, "^m68_p2\\[") | stringr::str_detect(par, "^m69_p1\\[") | stringr::str_detect(par, "^m69_p2\\[") | stringr::str_detect(par, "^m68_p3_2\\[") | stringr::str_detect(par, "^hpm_l_2_p3\\[")) %>%
                dplyr::mutate(
                  base_par =
                    dplyr::case_when(
                      stringr::str_detect(par, "^l_2\\[") ~ "l_2",
                      stringr::str_detect(par, "^k_2\\[") ~ "k_2",
                      stringr::str_detect(par, "^alpha_2\\[") ~ "alpha_2",
                      stringr::str_detect(par, "^m_rep\\[") ~ "m_rep",
                      stringr::str_detect(par, "^hpm_k_2\\[") ~ "hpm_k_2",
                      stringr::str_detect(par, "^hpm_l_2\\[") ~ "hpm_l_2",
                      stringr::str_detect(par, "^m68_p3_2\\[") ~ "m68_p3_2",
                      stringr::str_detect(par, "^hpm_l_2_p3\\[") ~ "hpm_l_2_p3",
                      TRUE ~ stringr::str_remove(par, pattern = "\\[\\d+\\]$")
                    )
                ) %>%
                dplyr::group_by(base_par) %>%
                dplyr::summarise(
                  id_fit = unique(id_fit),
                  dplyr::across(dplyr::starts_with("MCSE_"), function(.y) max(.y, na.rm = TRUE)),
                  .groups = "drop"
                )
            })
        )
    )

  res

}



#### RMSE ####

#' Computes the RMSE for one model
#'
#' @param x_stan_draws A data frame with omdel MCMC draws.
#'
#' @export
hpmd_get_rmse_one_model <- function(x_stan_draws, ...) {

  .dots <- list(...)

  res <-
    x_stan_draws %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::filter(index_hpm)

  if(.dots$has_cross_validation) {
    res <-
      res %>%
      dplyr::filter(ifelse(id_cross_validation_block_extraction == hpm_cross_validation_id_block, TRUE, FALSE))
  }

  res <-
    res %>%
    dplyr::mutate(
      rmse_total = sqrt(posterior::rvar_mean((k_2 - hpm_k_2_rep)^2, na.rm = TRUE)),
      rmse_for_cv_folds = sqrt(posterior::rvar_mean((k_2 - hpm_k_2_rep)[! is.na(hpm_cross_validation_id_block)]^2, na.rm = TRUE)),
      prediction_difference = k_2 - hpm_k_2_rep,
      prediction_difference_total = posterior::rvar_mean(prediction_difference, na.rm = TRUE)
    ) %>%
    dplyr::group_by(taxon_rank_value, taxon_rank_value_pretty) %>%
    dplyr::summarise(
      rmse_total = rmse_total[[1]],
      rmse_for_cv_folds = rmse_for_cv_folds[! is.na(rmse_for_cv_folds)][[1]],
      rmse = sqrt(posterior::rvar_mean((k_2 - hpm_k_2_rep)^2, na.rm = TRUE)),
      prediction_difference_total = prediction_difference_total[[1]],
      prediction_difference = posterior::rvar_mean(prediction_difference, na.rm = TRUE),
      k_2_sd = posterior::rvar_sd(k_2, na.rm = TRUE),
      k_2_mean = posterior::rvar_mean(k_2, na.rm = TRUE),
      n = length(k_2),
      n_id_citation = length(unique(id_citation)),
      list_id_citation = list(unique(id_citation)),
      .groups = "drop"
    )

  res

}



#' Computes the RMSE for one model
#'
#' @param x_stan_draws A data frame or list of data frames (if `has_cross_validation = TRUE`).
#'
#' @export
hpmd_get_rmse <- function(x_stan_draws, ...) {

  .dots <- list(...)

  res <-
    hpmd_get_rmse_one_model(
      x_stan_draws = x_stan_draws,
      ...
    )

  res

}
