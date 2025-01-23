#### Helper functions ####

#' Helper function to create simulation grids for simulations 3 and 4
#'
#' @param alpha_2_null A numeric value. The value for `alpha_2` for the assumed
#' null hypothesis.
#'
#' @param alpha_2_alternative A numeric value. The value for `alpha_2` for the
#' assumed alternative hypothesis.
#'
#' @export
hpmd_get_newdata_simulations_3_and_4 <- function(alpha_2_null = 2, alpha_2_alternative, x_stan_fit_4 = x_stan_fit_4, x_stan_data_4 = x_stan_data_4, hpmd_data = hpmd_data, sample_size = 100) {

  # sdm_mm31_1_parameters_standard_vs_estimated_p3_newdata
  newdata <-
    tibble::tibble(
      incubation_duration = c(seq(from = 0, to = 20, by = 1), seq(from = 21, to = 100, by = 5)),
      m0 = 1,
      layer_degree_of_saturation_1 = 0.6,
      layer_water_table_depth_to_surface_1 = 40,
      sample_depth_lower = 20,
      hpm_taxon_rank_value = "Sphagnum angustifolium"
    )

  newdata_1 <-
    dplyr::bind_rows(
      # S. angustifolium
      dplyr::bind_cols(
        newdata %>%
          dplyr::mutate(
            alpha_2 = alpha_2_null
          ) %>%
          hpmd_predict_fit_4(
            x_stan_fit_4 = x_stan_fit_4,
            x_stan_data_4 = x_stan_data_4,
            hpmd_data = hpmd_data
          ) %>%
          dplyr::rename(
            mass_relative_mass_mu_alpha_2 = "mass_relative_mass_mu"
          ),
        newdata %>%
          dplyr::mutate(
            alpha_2 = alpha_2_alternative
          ) %>%
          hpmd_predict_fit_4(
            x_stan_fit_4 = x_stan_fit_4,
            x_stan_data_4 = x_stan_data_4,
            hpmd_data = hpmd_data
          ) %>%
          dplyr::select(mass_relative_mass_mu) %>%
          dplyr::rename(
            mass_relative_mass_mu_alpha_1.5 = "mass_relative_mass_mu"
          )
      ),
      # S. fuscum
      dplyr::bind_cols(
        newdata %>%
          dplyr::mutate(
            hpm_taxon_rank_value = "Sphagnum fuscum",
            alpha_2 = alpha_2_null
          ) %>%
          hpmd_predict_fit_4(
            x_stan_fit_4 = x_stan_fit_4,
            x_stan_data_4 = x_stan_data_4,
            hpmd_data = hpmd_data
          ) %>%
          dplyr::rename(
            mass_relative_mass_mu_alpha_2 = "mass_relative_mass_mu"
          ),
        newdata %>%
          dplyr::mutate(
            hpm_taxon_rank_value = "Sphagnum fuscum",
            alpha_2 = alpha_2_alternative
          ) %>%
          hpmd_predict_fit_4(
            x_stan_fit_4 = x_stan_fit_4,
            x_stan_data_4 = x_stan_data_4,
            hpmd_data = hpmd_data
          ) %>%
          dplyr::select(mass_relative_mass_mu) %>%
          dplyr::rename(
            mass_relative_mass_mu_alpha_1.5 = "mass_relative_mass_mu"
          )
      )
    ) %>%
    dplyr::mutate(
      mass_relative_mass_mu_diff = mass_relative_mass_mu_alpha_2 - mass_relative_mass_mu_alpha_1.5
    )



  ## part 2
  # sdm_experiment_alpha_p2_newdata <-
  newdata <-
    tibble::tibble(
      incubation_duration = rep(5, sample_size),
      m0 = 1,
      layer_degree_of_saturation_1 = 0.6,
      layer_water_table_depth_to_surface_1 = 40,
      sample_depth_lower = 20,
      hpm_taxon_rank_value = "Sphagnum angustifolium"
    )

  new_column_name_alpha_2 <- c(mass_relative_mass_average_alpha_1.5 = paste0("mass_relative_mass_average_alpha_", alpha_2_alternative))

  newdata_2 <-
    dplyr::bind_rows(
      # S. angustifolium
      newdata %>%
        dplyr::mutate(
          alpha_2 = alpha_2_null
        ) %>%
        hpmd_predict_fit_4(
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data
        ) %>%
        dplyr::mutate(
          n = seq_len(nrow(.)),
          mass_relative_mass_average = cumsum(mass_relative_mass)/n
        ),
      newdata %>%
        dplyr::mutate(
          alpha_2 = alpha_2_alternative
        ) %>%
        hpmd_predict_fit_4(
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data
        ) %>%
        dplyr::mutate(
          n = seq_len(nrow(.)),
          mass_relative_mass_average = cumsum(mass_relative_mass)/n,
        ),
      # S. fuscum
      newdata %>%
        dplyr::mutate(
          alpha_2 = alpha_2_null,
          hpm_taxon_rank_value = "Sphagnum fuscum"
        ) %>%
        hpmd_predict_fit_4(
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data
        ) %>%
        dplyr::mutate(
          n = seq_len(nrow(.)),
          mass_relative_mass_average = cumsum(mass_relative_mass)/n
        ),
      newdata %>%
        dplyr::mutate(
          alpha_2 = alpha_2_alternative,
          hpm_taxon_rank_value = "Sphagnum fuscum"
        ) %>%
        hpmd_predict_fit_4(
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data
        ) %>%
        dplyr::mutate(
          n = seq_len(nrow(.)),
          mass_relative_mass_average = cumsum(mass_relative_mass)/n,
        )
    ) %>%
    dplyr::select(hpm_taxon_rank_value, incubation_duration, layer_degree_of_saturation_1, layer_water_table_depth_1, sample_depth_lower, mass_relative_mass_mu, mass_relative_mass_average, n, alpha_2) %>%
    tidyr::pivot_wider(
      names_from = alpha_2,
      values_from = c(mass_relative_mass_mu, mass_relative_mass_average),
      names_prefix = "alpha_"
    ) %>%
    dplyr::rename(
      dplyr::any_of(new_column_name_alpha_2)
    ) %>%
    dplyr::mutate(
      mass_relative_mass_diff = mass_relative_mass_average_alpha_2 - mass_relative_mass_average_alpha_1.5,
      pr_mass_relative_mass_diff_smaller_2 =
        if(alpha_2_null >= alpha_2_alternative) {
          posterior::Pr(mass_relative_mass_diff > 0)
        } else {
          posterior::Pr(mass_relative_mass_diff < 0)
        },
      index =
        {
          res <- rep(FALSE, nrow(.))

          index <- pr_mass_relative_mass_diff_smaller_2 >= 0.95 & hpm_taxon_rank_value == "Sphagnum fuscum"
          if(sum(index) > 0) {
            res[which(index)[[1]]] <- TRUE
          }
          index <- pr_mass_relative_mass_diff_smaller_2 >= 0.99 & hpm_taxon_rank_value == "Sphagnum fuscum"
          if(sum(index) > 0) {
            res[which(index)[[1]]] <- TRUE
          }
          index <- pr_mass_relative_mass_diff_smaller_2 >= 0.95 & hpm_taxon_rank_value == "Sphagnum angustifolium"
          if(sum(index) > 0) {
            res[which(index)[[1]]] <- TRUE
          }
          index <- pr_mass_relative_mass_diff_smaller_2 >= 0.99 & hpm_taxon_rank_value == "Sphagnum angustifolium"
          if(sum(index) > 0) {
            res[which(index)[[1]]] <- TRUE
          }
          res
        }
    )

  list(
    newdata_1 = newdata_1,
    newdata_2 = newdata_2
  )

}



#' Helper function to create plots for simulations 3 and 4
#'
#' @param newdata
#'
#' @export
hpmd_get_plot_for_simulations_3_and_4 <- function(newdata, file_plot = file_plot) {

  # sdm_experiment_alpha_p1
  p1 <-
    newdata$newdata_1 %>%
    dplyr::mutate(
      hpm_taxon_rank_value = paste0("<i>", hpm_taxon_rank_value, "</i>")
    ) %>%
    dplyr::select(mass_relative_mass_mu_diff, incubation_duration, hpm_taxon_rank_value) %>%
    ggplot(aes(ydist = mass_relative_mass_mu_diff * 100, x = incubation_duration)) +
    ggdist::stat_lineribbon() +
    geom_hline(yintercept = 0, color = "grey50") +
    scale_fill_brewer() +
    labs(
      y = "Difference of average remaining masses (%)",
      x = expression("Incubation duration (yr)")
    ) +
    facet_wrap(~ hpm_taxon_rank_value) +
    guides(
      fill = guide_legend(title = "Confidence interval")
    ) +
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.text = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      legend.position = "bottom"
    )

  # sdm_experiment_alpha_p2
  p2 <-
    newdata$newdata_2 %>%
    dplyr::mutate(
      hpm_taxon_rank_value = paste0("<i>", hpm_taxon_rank_value, "</i>")
    ) %>%
    dplyr::select(mass_relative_mass_diff, n, hpm_taxon_rank_value) %>%
    ggplot() +
    ggdist::stat_lineribbon(aes(ydist = mass_relative_mass_diff * 100, x = n)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_vline(
      data =
        newdata$newdata_2 %>%
        dplyr::mutate(
          hpm_taxon_rank_value = paste0("<i>", hpm_taxon_rank_value, "</i>")
        ) %>%
        dplyr::filter(index),
      mapping = aes(xintercept = n),
      color = "grey50"
    ) +
    scale_fill_brewer() +
    labs(
      y = "Difference of remaining masses (%)",
      x = expression("Sample size")
    ) +
    facet_wrap(~ hpm_taxon_rank_value) +
    guides(
      fill = guide_legend(title = "Confidence interval")
    ) +
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.text = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      legend.position = "bottom"
    )

  ### combine
  res_plot <-
    list(
      p1,
      p2
    ) %>%
    patchwork::wrap_plots(ncol = 1L) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    ) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(
    file_plot,
    plot = res_plot,
    width = 8.5, height = 7, dpi = 300,
    device = cairo_pdf
  )

}


#### Simulation grids ####


#' Creates template files for simulations
#'
#' @param id_simulation Integer value. The simulation for which to compute the
#' newdata data frame.
#'
#' @export
hpmd_get_newdata_for_simulations <- function(id_simulation, x_stan_fit_4 = hpmd_stan_fit_4, x_stan_data_4 = hpmd_data_stan_4, hpmd_data = hpmd_data) {

  # standard HPM parameter values
  hpmd_hpm_standard_parameter_values <-
    hpmd_get_hpm_standard_parameter_values()

  file_name <- paste0("_targets_rvars/hpmd_simulation_1_newdata_", id_simulation, ".rds")

  res <-
    switch(
      id_simulation,
      "1" = {
        # sdm_mm31_1_parameters_standard_vs_estimated_p1_newdata
        tibble::tibble(
          incubation_duration = 5,
          m0 = 1,
          layer_total_porosity_1 = 0.7,
          minimum_water_content_at_surface_1 = 0.05,
          layer_water_table_depth_to_surface_1 = 40,
          sample_depth_lower = c(seq(from = 0, to = 40, length.out = 11), seq(from = 40, to = 140, length.out = 11)[-1]),
          hpm_taxon_rank_value = "Sphagnum fuscum"
        ) %>%
          dplyr::mutate(
            layer_degree_of_saturation_1 =
              purrr::map_dbl(seq_along(sample_depth_lower), function(i) {
                m7(
                  layer_depth_midpoint_1 = sample_depth_lower[[i]],
                  layer_total_porosity_1 = layer_total_porosity_1[[i]],
                  water_table_depth_to_surface_1 = layer_water_table_depth_to_surface_1[[i]],
                  minimum_water_content_at_surface_1 = minimum_water_content_at_surface_1[[i]]
                ) / layer_total_porosity_1[[i]]
              })
          )
      },
      "2" = {
        # sdm_mm31_1_parameters_standard_vs_estimated_p2_newdata
        tibble::tibble(
          incubation_duration = 5,
          m0 = 1,
          layer_total_porosity_1 = 0.7,
          minimum_water_content_at_surface_1 = 0.05,
          layer_water_table_depth_to_surface_1 = 8,
          sample_depth_lower = c(seq(from = 0, to = 8, length.out = 11), seq(from = 8, to = 140, length.out = 11)[-1]),
          hpm_taxon_rank_value = "Sphagnum fallax"
        ) %>%
          dplyr::mutate(
            layer_degree_of_saturation_1 =
              purrr::map_dbl(seq_along(sample_depth_lower), function(i) {
                m7(
                  layer_depth_midpoint_1 = sample_depth_lower[[i]],
                  layer_total_porosity_1 = layer_total_porosity_1[[i]],
                  water_table_depth_to_surface_1 = layer_water_table_depth_to_surface_1[[i]],
                  minimum_water_content_at_surface_1 = minimum_water_content_at_surface_1[[i]]
                ) / layer_total_porosity_1[[i]]
              })
          )
      },
      "3" = {
        hpmd_get_newdata_simulations_3_and_4(
          alpha_2_null = 2,
          alpha_2_alternative = 1.5,
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data,
          sample_size = 100
        )
      },
      "4" = {
        hpmd_get_newdata_simulations_3_and_4(
          alpha_2_null = 2,
          alpha_2_alternative = 2.5,
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data,
          sample_size = 100
        )
      },
      "5" = {
        # same as simulation 2, but with standard value for m69_p1
        tibble::tibble(
          incubation_duration = 5,
          m0 = 1,
          layer_total_porosity_1 = 0.7,
          minimum_water_content_at_surface_1 = 0.05,
          layer_water_table_depth_to_surface_1 = 8,
          sample_depth_lower = c(seq(from = 0, to = 8, length.out = 11), seq(from = 8, to = 140, length.out = 11)[-1]),
          hpm_taxon_rank_value = "Sphagnum fallax",
          m69_p1 =
            hpmd_hpm_standard_parameter_values %>%
            dplyr::filter(variable == "m69_p1") %>%
            dplyr::pull(value)
        ) %>%
          dplyr::mutate(
            layer_degree_of_saturation_1 =
              purrr::map_dbl(seq_along(sample_depth_lower), function(i) {
                m7(
                  layer_depth_midpoint_1 = sample_depth_lower[[i]],
                  layer_total_porosity_1 = layer_total_porosity_1[[i]],
                  water_table_depth_to_surface_1 = layer_water_table_depth_to_surface_1[[i]],
                  minimum_water_content_at_surface_1 = minimum_water_content_at_surface_1[[i]]
                ) / layer_total_porosity_1[[i]]
              })
          )
      },
      "6" = {
        res <-
          tibble::tibble(
            incubation_duration = c(1:5, seq(10, 50, by = 5)),
            m0 = 1,
            layer_total_porosity_1 = 0.7,
            layer_degree_of_saturation_1 = 0.6,
            layer_water_table_depth_to_surface_1 = 30,
            sample_depth_lower = 20,
            hpm_taxon_rank_value = "Sphagnum fallax"
          )

        res <-
          dplyr::bind_rows(
            res,
            res %>%
              dplyr::mutate(
                sample_depth_lower = layer_water_table_depth_to_surface_1 + 20,
                layer_degree_of_saturation_1 = 1.0
              ),
            res %>%
              dplyr::mutate(
                hpm_taxon_rank_value = "Sphagnum fuscum"
              ),
            res %>%
              dplyr::mutate(
                sample_depth_lower = layer_water_table_depth_to_surface_1 + 20,
                layer_degree_of_saturation_1 = 1.0,
                hpm_taxon_rank_value = "Sphagnum fuscum"
              )
          )

        res

      }
    )

  saveRDS_rvars(res, file_name)

  file_name

}



#### Plots ####


#' Plots for the simulations
#'
#' @export
hpmd_get_plots_for_simulations <- function(id_simulation, newdata, x_stan_fit_4 = hpmd_stan_fit_4, x_stan_data_4 = hpmd_data_stan_4, hpmd_data = hpmd_data) {

  # standard HPM parameter values
  hpmd_hpm_standard_parameter_values <-
    hpmd_get_hpm_standard_parameter_values()

  file_plot <- paste0("figures/hpmd_simulation_1_plot_1_", id_simulation, ".pdf")

  switch(
    id_simulation,
    "1" = {

      # sdm_mm31_1_parameters_standard_vs_estimated_p1
      res_plot <-
        dplyr::bind_rows(
          # Wopt (m69_p1)
          dplyr::bind_cols(
            newdata %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ),
            newdata %>%
              dplyr::mutate(
                m69_p1 = posterior::as_rvar(hpmd_hpm_standard_parameter_values$value[hpmd_hpm_standard_parameter_values$variable == "m69_p1"]),
                # m69_p2 = posterior::as_rvar(2.31),
                # m68_p1 = 0.001,
                # m68_p2 = 0.3
              ) %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ) %>%
              dplyr::select(hpm_k_2) %>%
              dplyr::rename(
                hpm_k_2_std = "hpm_k_2"
              )
          ) %>%
            dplyr::mutate(
              parameter = "m69_p1",
              parameter_pretty = hpmd_hpm_standard_parameter_values$variable_pretty[hpmd_hpm_standard_parameter_values$variable == "m69_p1"]
            ),
          # c1 (m69_p2)
          dplyr::bind_cols(
            newdata %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ),
            newdata %>%
              dplyr::mutate(
                # m69_p1 = posterior::as_rvar(0.45),
                m69_p2 = posterior::as_rvar(hpmd_hpm_standard_parameter_values$value[hpmd_hpm_standard_parameter_values$variable == "m69_p2"]),
                # m68_p1 = 0.001,
                # m68_p2 = 0.3
              ) %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ) %>%
              dplyr::select(hpm_k_2) %>%
              dplyr::rename(
                hpm_k_2_std = "hpm_k_2"
              )
          ) %>%
            dplyr::mutate(
              parameter = "m69_p2",
              parameter_pretty = hpmd_hpm_standard_parameter_values$variable_pretty[hpmd_hpm_standard_parameter_values$variable == "m69_p2"]
            ),
          # fmin (m68_p1)
          dplyr::bind_cols(
            newdata %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ),
            newdata %>%
              dplyr::mutate(
                # m69_p1 = posterior::as_rvar(0.45),
                # m69_p2 = posterior::as_rvar(2.31),
                m68_p1 = posterior::as_rvar(hpmd_hpm_standard_parameter_values$value[hpmd_hpm_standard_parameter_values$variable == "m68_p1"]),
                # m68_p2 = 0.3
              ) %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ) %>%
              dplyr::select(hpm_k_2) %>%
              dplyr::rename(
                hpm_k_2_std = "hpm_k_2"
              )
          ) %>%
            dplyr::mutate(
              parameter = "m68_p1",
              parameter_pretty = hpmd_hpm_standard_parameter_values$variable_pretty[hpmd_hpm_standard_parameter_values$variable == "m68_p1"]
            ),
          # c2 (m68_p2)
          dplyr::bind_cols(
            newdata %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ),
            newdata %>%
              dplyr::mutate(
                # m69_p1 = posterior::as_rvar(0.45),
                # m69_p2 = posterior::as_rvar(2.31),
                # m68_p1 = 0.001,
                m68_p2 = posterior::as_rvar(hpmd_hpm_standard_parameter_values$value[hpmd_hpm_standard_parameter_values$variable == "m68_p2"]),
              ) %>%
              hpmd_predict_fit_4(
                x_stan_fit_4 = x_stan_fit_4,
                x_stan_data_4 = x_stan_data_4,
                hpmd_data = hpmd_data
              ) %>%
              dplyr::select(hpm_k_2) %>%
              dplyr::rename(
                hpm_k_2_std = "hpm_k_2"
              )
          ) %>%
            dplyr::mutate(
              parameter = "m68_p2",
              parameter_pretty = hpmd_hpm_standard_parameter_values$variable_pretty[hpmd_hpm_standard_parameter_values$variable == "m68_p2"]
            )
        ) %>%
        dplyr::mutate(
          x = layer_water_table_depth_1,
          parameter_pretty =
            parameter_pretty %>%
            factor(levels = unique(.))
        ) %>%
        dplyr::select(x, hpm_k_2_std, hpm_k_2, parameter_pretty) %>%
        ggplot(aes(ydist = hpm_k_2 - hpm_k_2_std, x = x)) +
        ggdist::stat_lineribbon() +
        geom_hline(yintercept = 0, color = "grey50") +
        scale_fill_brewer() +
        labs(
          y = paste0("<i>k<sub>0,modified</sub></i>(", hpmd_model_id_to_name_2("3-4"), ") - <i>k<sub>0,standard</sub></i>(", hpmd_model_id_to_name_2("3-4"), ") (yr<sup>-1</sup>)"),
          x = expression("Water table depth below litter (cm)")
        ) +
        facet_wrap(~ parameter_pretty, ncol = 2L, scales = "free_x") +
        guides(
          fill = guide_legend(title = "Confidence interval")
        ) +
        theme(
          axis.title.y = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(size = 12),
          strip.background.x = element_blank(),
          legend.position = "bottom"
        )

      # export
      ggsave(
        file_plot,
        plot = res_plot,
        width = 5.5, height = 5, dpi = 300,
        device = cairo_pdf
      )

    },
    "2" = {

      # sdm_mm31_1_parameters_standard_vs_estimated_p2
      res <- readRDS_rvars(tar_read(hpmd_simulation_1_newdata_1))
      res_plot <-
        dplyr::bind_rows(
          res %>%
            dplyr::mutate(
              microhabitat = "Hummock",
              wopt = paste0("<i>W</i><sub><i>opt</i>, ", hpmd_model_id_to_name_2("3-4"), "</sub>")
            ) %>%
            hpmd_predict_fit_4(
              x_stan_fit_4 = x_stan_fit_4,
              x_stan_data_4 = x_stan_data_4,
              hpmd_data = hpmd_data
            ),
          res %>%
            dplyr::mutate(
              microhabitat = "Hummock",
              m69_p1 = posterior::as_rvar(hpmd_hpm_standard_parameter_values$value[hpmd_hpm_standard_parameter_values$variable == "m69_p1"]),
              wopt = paste0("<i>W</i><sub><i>opt</i>, standard</sub></i>")
            ) %>%
            hpmd_predict_fit_4(
              x_stan_fit_4 = x_stan_fit_4,
              x_stan_data_4 = x_stan_data_4,
              hpmd_data = hpmd_data
            )
        ) %>%
        dplyr::mutate(
          wopt = factor(wopt, levels = c(paste0("<i>W</i><sub><i>opt</i>, standard</sub></i>"), paste0("<i>W</i><sub><i>opt</i>, ", hpmd_model_id_to_name_2("3-4"), "</sub>")))
        ) %>%
        dplyr::filter(sample_depth_lower <= 50) %>%
        dplyr::select(sample_depth_lower, layer_water_table_depth_to_surface_1, hpm_k_2, wopt) %>%
        ggplot(aes(ydist = hpm_k_2, x = sample_depth_lower)) +
        ggdist::stat_lineribbon() +
        geom_vline(aes(xintercept = layer_water_table_depth_to_surface_1), color = "grey50") +
        geom_hline(aes(yintercept = 0.0), color = "grey50") +
        scale_fill_brewer() +
        facet_wrap(~ wopt, scales = "free_x") +
        coord_flip() +
        scale_x_reverse() +
        labs(
          y = paste0("<i>k<sub>0</sub></i><sub>, modified</sub>(", hpmd_model_id_to_name_2("3-4"), ") (yr<sup>-1</sup>)"),
          x = expression("Sample depth (cm)")
        ) +
        guides(
          fill = guide_legend(title = "Confidence interval")
        ) +
        theme(
          axis.title.y = ggtext::element_markdown(),
          axis.title.x = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(size = 12),
          strip.background.x = element_blank(),
          legend.position = "bottom"
        )

      # export
      ggsave(
        file_plot,
        plot = res_plot,
        width = 3.5, height = 4, dpi = 300,
        device = cairo_pdf
      )

    },
    "3" = {
      hpmd_get_plot_for_simulations_3_and_4(newdata = newdata, file_plot = file_plot)
    },
    "4" = {
      hpmd_get_plot_for_simulations_3_and_4(newdata = newdata, file_plot = file_plot)
    },
    "5" = {

      res_plot <-
        dplyr::bind_rows(
          newdata %>%
            dplyr::mutate(
              microhabitat = "Hollow",
              wopt = paste0("<i>W</i><sub><i>opt</i>, ", hpmd_model_id_to_name_2("3-4"), "</sub>")
            ) %>%
            hpmd_predict_fit_4(
              x_stan_fit_4 = x_stan_fit_4,
              x_stan_data_4 = x_stan_data_4,
              hpmd_data = hpmd_data
            ),
          newdata %>%
            dplyr::mutate(
              m69_p1 =
                hpmd_hpm_standard_parameter_values %>%
                dplyr::filter(variable == "m69_p1") %>%
                dplyr::pull(value),
              microhabitat = "Hollow",
              wopt = paste0("<i>W</i><sub><i>opt</i>, standard</sub></i>")
            ) %>%
            hpmd_predict_fit_4(
              x_stan_fit_4 = x_stan_fit_4,
              x_stan_data_4 = x_stan_data_4,
              hpmd_data = hpmd_data
            )
        )  %>%
        dplyr::mutate(
          wopt = factor(wopt, levels = c(paste0("<i>W</i><sub><i>opt</i>, standard</sub></i>"), paste0("<i>W</i><sub><i>opt</i>, ", hpmd_model_id_to_name_2("3-4"), "</sub>")))
        ) %>%
        dplyr::filter(sample_depth_lower <= 50) %>%
        dplyr::select(sample_depth_lower, layer_water_table_depth_to_surface_1, hpm_k_2, wopt) %>%
        ggplot(aes(ydist = hpm_k_2, x = sample_depth_lower)) +
        ggdist::stat_lineribbon() +
        geom_vline(aes(xintercept = layer_water_table_depth_to_surface_1), color = "grey50") +
        scale_fill_brewer() +
        facet_wrap(~ wopt, scales = "free_x") +
        coord_flip() +
        scale_x_reverse() +
        labs(
          y = paste0("<i>k<sub>0</sub></i><sub>, modified</sub>(", hpmd_model_id_to_name_2("3-4"), ") (yr<sup>-1</sup>)"),
          x = expression("Sample depth (cm)")
        ) +
        guides(
          fill = guide_legend(title = "Confidence interval")
        ) +
        theme(
          axis.title.y = ggtext::element_markdown(),
          axis.title.x = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(size = 12),
          strip.background.x = element_blank(),
          legend.position = "bottom"
        )

      # export
      ggsave(
        file_plot,
        plot = res_plot,
        width = 3.5, height = 4, dpi = 300,
        device = cairo_pdf
      )

    },
    "6" = {
      res <-
        newdata %>%
        hpmd_predict_fit_4(
          x_stan_fit_4 = x_stan_fit_4,
          x_stan_data_4 = x_stan_data_4,
          hpmd_data = hpmd_data
        )

      # average remaining mass
      p1 <-
        res %>%
        ggplot(aes(ydist = mass_relative_mass_mu * 100, x = incubation_duration)) +
        ggdist::stat_lineribbon(aes(fill = hpm_taxon_rank_value), alpha = 0.2, .width = c(0.5, 0.8, 0.95)) +
        facet_wrap(~ paste0(layer_water_table_depth_1, " cm")) +
        labs(
          y = "Fraction of initial mass (%)",
          x = "Incubation duration (yr)",
          title = "Average remaining mass"
        ) +
        guides(fill = guide_legend(title = "Species")) +
        scale_fill_manual(values = c("salmon", "steelblue")) +
        theme(
          strip.background.x = element_blank()
        )

      # estimate for one sample
      p2 <-
        res %>%
        ggplot(aes(ydist = mass_relative_mass * 100, x = incubation_duration)) +
        ggdist::stat_lineribbon(aes(fill = hpm_taxon_rank_value), alpha = 0.2, .width = c(0.5, 0.8, 0.95)) +
        facet_wrap(~ paste0(layer_water_table_depth_1, " cm")) +
        labs(
          y = "Fraction of initial mass (%)",
          x = "Incubation duration (yr)",
          title = "Remaining mass for one sample"
        ) +
        guides(fill = guide_legend(title = "Species")) +
        scale_fill_manual(values = c("salmon", "steelblue")) +
        theme(
          strip.background.x = element_blank()
        )

      ## combine
      res_plot <-
        list(
          p1,
          p2
        ) %>%
        patchwork::wrap_plots(ncol = 2L) +
        patchwork::plot_annotation(
          tag_levels = c('a'),
          tag_prefix = '(',
          tag_suffix = ')'
        ) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

      ggsave(
        file_plot,
        plot = res_plot,
        width = 8.5, height = 4, dpi = 300,
        device = cairo_pdf
      )


    }
  )

  file_plot

}
