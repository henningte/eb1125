#### Plots ####

#' k$_0$ estimated from litterbag data versus k$_0$ predicted by the HPM
#'
#' @param x_stan_draws A list of character values poiting to files storing data
#' frames with MCMC draws from the models.
#'
#' @export
hpmd_get_plot_1 <- function(x_stan_draws) {

  # combine data
  res <-
    x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars) %>%
    dplyr::left_join(
      hpmd_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model =
        model_name |>
        hpmd_model_id_to_name() |>
        factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name()),
      hpm_microhabitat =
        factor(hpm_microhabitat, levels = c("Hollow", "Lawn", "Hummock"))
    )

  res_plot <-
    res %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    ggplot() +
    ggdist::stat_interval(aes(ydist = k_2, x = mean(hpm_k_2_rep)), .width = 0.95, interval_size = 0.5, interval_colour = "grey", show.legend = FALSE) +
    ggdist::stat_interval(aes(xdist = hpm_k_2_rep, y = mean(k_2)), .width = 0.95, interval_size = 0.5, interval_colour = "grey", show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    geom_point(aes(x = mean(hpm_k_2_rep), y = mean(k_2), fill = hpm_microhabitat), shape = 21, size = 2) +
    labs(
      y = expression("Estimated from litterbag data (yr"^{-1}*")"),
      x = expression("Predicted by HPM (modifications) (yr"^{-1}*")")
    ) +
    scale_fill_manual(values = c("steelblue", "yellowgreen", "salmon")) +
    facet_wrap(~ id_model, nrow = 1L) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12)
    ) +
    guides(
      fill = guide_legend(title = "Microhabitat", override.aes = list(size = 3))
    ) +
    coord_cartesian(ylim = c(NA, 0.5))

  file_plot <- "figures/hpmd_plot_1.pdf"

  # export
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 3.5, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}



#' k$_0$ estimated from litterbag data or predicted by the HPM versus water table depth
#'
#' @param x_stan_draws
#'
#' @export
hpmd_get_plot_2 <- function(x_stan_draws) {

  # function to add custom breaks. Modified from https://stackoverflow.com/a/55980394
  equal_breaks <- function(n = 3, s = 0.05, ..., digits){
    function(x){
      # rescaling
      d <- s * diff(range(x)) / (1+2*s)
      round(seq(min(x)+d, max(x)-d, length=n), digits = digits)
    }
  }

  # combine data
  res <-
    x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars) %>%
    dplyr::left_join(
      hpmd_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model =
        model_name |>
        hpmd_model_id_to_name() |>
        factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
    )

  res_plot <-
    res %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    dplyr::filter(
      taxon_rank_value %in% c("Sphagnum", "Sphagnum angustifolium", "Sphagnum balticum", "Sphagnum cuspidatum", "Sphagnum fuscum", "Sphagnum magellanicum")
    ) %>%
    dplyr::mutate(
      x =
        dplyr::case_when(
          id_model == "HPMf" ~ layer_water_table_depth_1,
          TRUE ~ layer_water_table_depth_to_surface_1 - sample_depth_upper
        )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("k_2", "hpm_k_2_rep")),
      names_to = "variable",
      values_to = "y"
    ) %>%
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "k_2" ~ "No",
          TRUE ~ "Yes"
        )
      ,
      y_mean = mean(y),
      y_lower = quantile(y, probs = 0.025)[1, ],
      y_upper = quantile(y, probs = 0.975)[1, ],
      x_mean = mean(x),
      x_lower = quantile(x, probs = 0.025)[1, ],
      x_upper = quantile(x, probs = 0.975)[1, ]
    ) %>%
    ggplot(aes(y = y_mean, x = x_mean)) +
    geom_errorbarh(aes(xmin = x_lower, xmax = x_upper), height = 0, color = "grey") +
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), width = 0, color = "grey") +
    geom_smooth(aes(color = variable, fill = variable), formula = y ~ x, method = "lm", se = FALSE) +
    geom_point(aes(fill = variable), shape = 21) +
    facet_grid(taxon_rank_value_pretty ~ id_model, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA), breaks = equal_breaks(n = 4, s = 0.05, digits = 2)) +
    scale_fill_manual(values = c("black", "grey")) +
    scale_color_manual(values = c("black", "grey")) +
    labs(
      y = expression("Decomposition rate (yr"^{-1}*")"),
      x = "Water table depth below litterbag (cm)"
    ) +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      strip.text.x = ggtext::element_markdown(size = 12),
      strip.text.y = ggtext::element_markdown(angle = 0, hjust = 0, size = 11),
      panel.spacing.y = unit(0.8, "lines")
    ) +
    guides(
      color = guide_legend(title = "Predicted with HPM?", override.aes = list(size = 3)),
      fill = guide_legend(title = "Predicted with HPM?")
    )

  file_plot <- "figures/hpmd_plot_2.pdf"

  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 7, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}


#' k$_0$ estimated by the HPM versus k$_0$ estimated by model 1-4 from the leaching project and l$_0$ estimated by the HPM versus l$_0$ estimated by model 1-4 from the leaching project
#'
#' @export
hpmd_get_plot_3 <- function(x_stan_draws) {

  # combine data
  res <-
    x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars)

  res <-
    res %>%
    dplyr::left_join(
      hpmd_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model =
        model_name |>
        hpmd_model_id_to_name() |>
        factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
    ) %>%
    dplyr::left_join(
      res %>%
        dplyr::filter(id_fit == 1) %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
        dplyr::select(id_sample_incubation_start, k_2, l_2) %>%
        dplyr::rename(
          k_2_id_fit_1 = "k_2",
          l_2_id_fit_1 = "l_2"
        ),
      by = "id_sample_incubation_start"
    )

  # initial leaching losses
  p1 <-
    res %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    dplyr::filter(! id_citation == "Hagemann.2015") %>%
    dplyr::mutate(
      l_2 = l_2 * 100,
      l_2_id_fit_1 = l_2_id_fit_1 * 100,
      y_mean = mean(l_2),
      y_lower = quantile(l_2, probs = 0.025)[1, ],
      y_upper = quantile(l_2, probs = 0.975)[1, ],
      x_mean = mean(l_2_id_fit_1),
      x_lower = quantile(l_2_id_fit_1, probs = 0.025)[1, ],
      x_upper = quantile(l_2_id_fit_1, probs = 0.975)[1, ],
      index_hpm =
        dplyr::case_when(
          index_hpm ~ "Yes",
          TRUE ~ "No"
        ) %>%
        factor(levels = c("Yes", "No"))
    ) %>%
    dplyr::arrange(dplyr::desc(index_hpm)) %>%
    ggplot(aes(y = y_mean, x = x_mean)) +
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), color = "grey", width = 0) +
    geom_errorbarh(aes(xmin = x_lower, xmax = x_upper), color = "grey", height = 0) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    geom_point(aes(fill = index_hpm), size = 2, shape = 21) +
    labs(
      title = expression("Initial leaching loss (mass-%)"),
      x = expression("Estimated with litterbag model (mass-%)"),
      y = expression("Estimated with HPM modifications (mass-%)")
    ) +
    scale_fill_manual(values = c("white", "black")) +
    facet_wrap(~ id_model) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12)
    ) +
    guides(
      fill = guide_legend(title = "Was modelled with HPM?", override.aes = list(size = 3))
    )

  # decomposition rates
  p2 <-
    res %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    dplyr::filter(! id_citation == "Hagemann.2015") %>%
    dplyr::mutate(
      y_mean = mean(k_2),
      y_lower = quantile(k_2, probs = 0.025)[1, ],
      y_upper = quantile(k_2, probs = 0.975)[1, ],
      x_mean = mean(k_2_id_fit_1),
      x_lower = quantile(k_2_id_fit_1, probs = 0.025)[1, ],
      x_upper = quantile(k_2_id_fit_1, probs = 0.975)[1, ],
      index_hpm =
        dplyr::case_when(
          index_hpm ~ "Yes",
          TRUE ~ "No"
        ) %>%
        factor(levels = c("Yes", "No"))
    ) %>%
    dplyr::arrange(dplyr::desc(index_hpm)) %>%
    ggplot(aes(y = y_mean, x = x_mean)) +
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), color = "grey", width = 0) +
    geom_errorbarh(aes(xmin = x_lower, xmax = x_upper), color = "grey", height = 0) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    geom_point(aes(fill = index_hpm), size = 2, shape = 21) +
    labs(
      title = expression("Decomposition rate (yr"^{-1}*")"),
      x = expression("Estimated with litterbag model (yr"^{-1}*")"),
      y = expression("Estimated with HPM modifications (yr"^{-1}*")")
    ) +
    scale_fill_manual(values = c("white", "black")) +
    facet_wrap(~ id_model) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12)
    ) +
    guides(
      fill = guide_legend(title = "Was modelled with HPM?", override.aes = list(size = 3))
    )

  res_plot <-
    list(p1, p2) %>%
    patchwork::wrap_plots(nrow = 1L, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    ) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # export
  file_plot <- "figures/hpmd_plot_3.pdf"
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 5, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}



#' HPM parameter estimates
#'
#' @param x_stan_draws
#'
#' @export
hpmd_get_plot_4 <- function(x_stan_draws) {

  # get parameter values
  hpmd_hpm_standard_parameter_values <-
    hpmd_get_hpm_standard_parameter_values()

  ### Are standard parameter values of the HPM compatible with our estimates?
  p1 <-
    x_stan_draws %>%
    dplyr::filter(! is.na(m69_p1)) %>%
    dplyr::slice(1) %>%
    dplyr::select(m69_p1, m69_p2, m68_p1, m68_p2) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        function(.x) list(attr(.x, "draws")[, 1])
      )
    ) %>%
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::left_join(
      hpmd_hpm_standard_parameter_values %>%
        dplyr::rename(
          value_default = "value"
        ) %>%
        dplyr::select(variable, value_default, variable_pretty),
      by = "variable"
    ) %>%
    tidyr::unnest(
      value
    ) %>%
    ggplot(aes(x = value)) +
    geom_density(fill = "grey") +
    geom_vline(aes(xintercept = value_default)) +
    labs(
      y = "Density",
      x = "Parameter value"
    ) +
    facet_wrap(~ variable_pretty, scales = "free", strip.position = "top") +
    theme(
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12),
      strip.placement = "outside"
    )

  ### Initial decomposition rates estimated for each species versus default values

  # labels
  p2_labels <-
    x_stan_draws %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(taxon_rank_value)) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::select(taxon_rank_value, taxon_rank_value_pretty, hpm_microhabitat, m68_p3_2) %>%
    dplyr::mutate(
      x_end = mean(m68_p3_2),
      x_start =
        dplyr::case_when(
          taxon_rank_value == "Sphagnum papillosum" ~ x_end + 0.035,
          taxon_rank_value == "Sphagnum balticum" ~ x_end - 0.05,
          taxon_rank_value == "Sphagnum fallax" ~ x_end + 0.07,
          taxon_rank_value == "Sphagnum fuscum" ~ x_end,
          taxon_rank_value == "Sphagnum" ~ x_end + 0.04,
          taxon_rank_value == "Sphagnum magellanicum" ~ x_end + 0.05,
          taxon_rank_value == "Sphagnum rubellum" ~ x_end + 0.09,
          taxon_rank_value == "Sphagnum angustifolium" ~ x_end + 0.17,
          taxon_rank_value == "Sphagnum russowii" ~ x_end + 0.12,
          taxon_rank_value == "Sphagnum majus" ~ x_end + 0.05,
          TRUE ~ x_end
        )
    )

  # plot
  p2 <-
    x_stan_draws %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(taxon_rank_value)) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::select(taxon_rank_value, taxon_rank_value_pretty, hpm_microhabitat, m68_p3_2) %>%
    dplyr::left_join(
      hpmd_hpm_standard_parameter_values %>%
        dplyr::rename(
          k0 = "value"
        ) %>%
        dplyr::select(hpm_microhabitat, k0),
      by = "hpm_microhabitat"
    ) %>%
    dplyr::mutate(
      m68_p3_2_mean = mean(m68_p3_2),
      m68_p3_2 =
        purrr::map(seq_along(m68_p3_2), function(i) tibble::tibble(m68_p3_2 = attr(m68_p3_2[[i]], "draws")[, 1]))
    ) %>%
    tidyr::unnest(
      m68_p3_2
    ) %>%
    dplyr::mutate(
      hpm_microhabitat = factor(hpm_microhabitat, levels = c("Hollow", "Lawn", "Hummock"))
    ) %>%
    ggplot() +
    geom_density(aes(x = m68_p3_2, group = taxon_rank_value), fill = "grey", alpha = 0.3) +
    geom_segment(aes(x = k0, xend = k0, y = 0, yend = 50)) +
    geom_segment(
      data = p2_labels,
      aes(x = x_start, xend = x_end, y = -13, yend = -2),
      color = "grey50"
    ) +
    ggtext::geom_richtext(
      data = p2_labels,
      aes(x = x_start, y = -15, label = taxon_rank_value_pretty),
      angle = 90, hjust = 1,
      fill = NA, label.color = NA, # remove background and outline
      label.padding = grid::unit(rep(0, 4), "pt") # remove padding
    ) +
    scale_y_continuous(limits = c(-90, NA)) +
    labs(
      y = "Density",
      x = "<i>k<sub>0, i</sub></i> (yr<sup>-1</sup>)"
    ) +
    facet_wrap(~ hpm_microhabitat) +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      strip.text.x = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown()
    )

  # combine
  res_plot <-
    list(p2, p1) %>%
    patchwork::wrap_plots(nrow = 1L, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    )

  # export
  file_plot <- paste0("figures/hpmd_plot_4_", x_stan_draws$id_fit[[1]], ".pdf")
  ggsave(
    file_plot,
    plot = res_plot,
    width = 11, height = 4.5, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}


#' HPM parameter values in cross-validation fits
#'
#' @export
hpmd_get_plot_5 <- function(x_stan_draws) {

  hpmd_hpm_standard_parameter_values <-
    hpmd_get_hpm_standard_parameter_values()

  ### HPM parameters except k_{0,i}
  p1 <-
    x_stan_draws %>%
    dplyr::filter(index_hpm & ifelse(hpm_cross_validation_id_block == id_cross_validation_block_extraction, TRUE, FALSE)) %>%
    dplyr::filter(! duplicated(paste0(id_cross_validation_block_extraction))) %>%
    dplyr::select(hpm_cross_validation_id_block, m68_p1, m68_p2, m69_p1, m69_p2) %>%
    tidyr::pivot_longer(
      ! dplyr::all_of("hpm_cross_validation_id_block"),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    dplyr::left_join(
      hpmd_hpm_standard_parameter_values %>%
        dplyr::select(variable, variable_pretty) %>%
        dplyr::rename(
          parameter = "variable",
          parameter_pretty = "variable_pretty"
        ),
      by = "parameter"
    ) %>%
    dplyr::mutate(
      hpm_cross_validation_id_block = as.factor(hpm_cross_validation_id_block)
    ) %>%
    ggplot(aes(xdist = value, group = hpm_cross_validation_id_block)) +
    ggdist::stat_slab(normalize = "panels", fill = NA, color = "black") +
    facet_wrap(~ parameter_pretty, scales = "free_x") +
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.text = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      y = "Density",
      x = "Parameter value"
    )

  ### k{0,i}

  ## for this, we need to extract new paramater values because d_mm32_1_yhat contains only the parameter values for the tested CV fold
  # extract predictions
  p2 <-
    x_stan_draws %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(paste0(hpm_microhabitat2, "_", id_cross_validation_block_extraction))) %>%
    dplyr::filter(taxon_rank_value %in% taxon_rank_value[! is.na(hpm_cross_validation_id_block)]) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    ggplot(aes(xdist = m68_p3_2, group = hpm_cross_validation_id_block)) +
    ggdist::stat_slab(normalize = "panels", fill = NA, color = "black") +
    facet_wrap(~ taxon_rank_value_pretty, scales = "free_x") +
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.text = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      y = "Density (-)",
      x = "<i>k<sub>0,i</sub></i> (yr<sup>-1</sup>)"
    )


  ## combine
  res_plot <-
    list(
      p2,
      p1
    ) %>%
    patchwork::wrap_plots(nrow = 1L, byrow = TRUE, widths = c(1, 0.7)) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    )

  ## export
  file_plot <- "figures/hpmd_plot_5.pdf"
  ggsave(
    file_plot,
    plot = res_plot,
    width = 10, height = 4.5, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}


#' $k_{0,i}$ for model 3-3 and model 3-4
#'
#' @export
hpmd_get_plot_6 <- function(x_stan_draws, hpmd_d_models) {

  res_plot <-
    purrr::map(x_stan_draws, function(.x) {
      .x %>%
        dplyr::filter(index_hpm) %>%
        dplyr::filter(! duplicated(hpm_microhabitat2)) %>%
        leaching_add_taxon_rank_value_pretty() %>%
        dplyr::mutate(
          id_model =
            hpmd_d_models$model_name[hpmd_d_models$id_fit == id_fit[[1]]] |>
            hpmd_model_id_to_name() |>
            factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
        ) %>%
        dplyr::select(hpm_microhabitat2, taxon_rank_value_pretty, m68_p3_2, id_model)
    }) %>%
    dplyr::bind_rows() %>%
    ggplot(aes(x = taxon_rank_value_pretty, color = id_model, ydist = m68_p3_2)) +
    ggdist::stat_pointinterval(.width = c(0.95)) +
    coord_flip() +
    scale_color_manual(values = c("black", "grey")) +
    labs(
      y = "Species",
      x = "<i>k</i><sub>0,i</sub> (yr<sup>-1</sup>)"
    ) +
    guides(
      color = guide_legend(title = "Model")
    ) +
    theme(
      axis.text.y = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.position = "bottom"
    )

  file_plot <- "figures/hpmd_plot_6.pdf"
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}



#' $l_0$ estimated from litterbag data or predicted by the HPM versus water table depth
#'
#' @param x_stan_draws
#'
#' @export
hpmd_get_plot_7 <- function(x_stan_draws) {

  # function to add custom breaks. Modified from https://stackoverflow.com/a/55980394
  equal_breaks <- function(n = 3, s = 0.05, ..., digits){
    function(x){
      # rescaling
      d <- s * diff(range(x)) / (1+2*s)
      round(seq(min(x)+d, max(x)-d, length=n), digits = digits)
    }
  }

  # combine data
  res <-
    x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars) %>%
    dplyr::left_join(
      hpmd_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model =
        model_name |>
        hpmd_model_id_to_name() |>
        factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
    )

  res_plot <-
    res %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    dplyr::filter(
      taxon_rank_value %in% c("Sphagnum", "Sphagnum angustifolium", "Sphagnum balticum", "Sphagnum cuspidatum", "Sphagnum fuscum", "Sphagnum magellanicum")
    ) %>%
    dplyr::mutate(
      x =
        dplyr::case_when(
          id_model == "HPMf" ~ layer_water_table_depth_1,
          TRUE ~ layer_water_table_depth_to_surface_1 - sample_depth_upper
        )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("l_2", "hpm_l_2_rep")),
      names_to = "variable",
      values_to = "y"
    ) %>%
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "l_2" ~ "No",
          TRUE ~ "Yes"
        ),
      y = y * 100,
      y_mean = mean(y),
      y_lower = quantile(y, probs = 0.025)[1, ],
      y_upper = quantile(y, probs = 0.975)[1, ],
      x_mean = mean(x),
      x_lower = quantile(x, probs = 0.025)[1, ],
      x_upper = quantile(x, probs = 0.975)[1, ]
    ) %>%
    ggplot(aes(y = y_mean, x = x_mean)) +
    geom_errorbarh(aes(xmin = x_lower, xmax = x_upper), height = 0, color = "grey") +
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), width = 0, color = "grey") +
    geom_smooth(aes(color = variable, fill = variable), formula = y ~ x, method = "lm", se = FALSE) +
    geom_point(aes(fill = variable), shape = 21) +
    facet_grid(taxon_rank_value_pretty ~ id_model, scales = "free_y") +
    # scale_y_continuous(limits = c(0, NA), breaks = equal_breaks(n = 4, s = 5, digits = 0)) +
    scale_fill_manual(values = c("black", "grey")) +
    scale_color_manual(values = c("black", "grey")) +
    labs(
      y = expression("Initial leaching loss (mass-%)"),
      x = "Water table depth below litterbag (cm)"
    ) +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      strip.text.x = ggtext::element_markdown(size = 12),
      strip.text.y = ggtext::element_markdown(angle = 0, hjust = 0, size = 11),
      panel.spacing.y = unit(0.8, "lines")
    ) +
    guides(
      color = guide_legend(title = "Predicted with HPM?", override.aes = list(size = 3)),
      fill = guide_legend(title = "Predicted with HPM?")
    )

  file_plot <- "figures/hpmd_plot_7.pdf"

  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 7, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}



#' k$_0$ estimated from litterbag data or predicted by the HPM versus water table depth for selected models and grouped by `id_citation`
#'
#' @param x_stan_draws A data frame with MCMC draws from one model.
#'
#' @export
hpmd_get_plot_8 <- function(x_stan_draws, file_plot) {

  # function to add custom breaks. Modified from https://stackoverflow.com/a/55980394
  equal_breaks <- function(n = 3, s = 0.05, ..., digits){
    function(x){
      # rescaling
      d <- s * diff(range(x)) / (1+2*s)
      round(seq(min(x)+d, max(x)-d, length=n), digits = digits)
    }
  }

  # combine data
  res <-
    x_stan_draws %>%
    dplyr::left_join(
      hpmd_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model =
        model_name |>
        hpmd_model_id_to_name() |>
        factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
    ) %>%
    dplyr::bind_rows(
      tibble::tibble(layer_water_table_depth_1 = numeric(), layer_water_table_depth_to_surface_1 = numeric())
    )

  res_plot <-
    res %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    dplyr::filter(
      taxon_rank_value %in% c("Sphagnum", "Sphagnum angustifolium", "Sphagnum balticum", "Sphagnum cuspidatum", "Sphagnum fuscum", "Sphagnum magellanicum")
    ) %>%
    dplyr::mutate(
      x =
        dplyr::case_when(
          id_model == "HPMf" ~ layer_water_table_depth_1,
          TRUE ~ layer_water_table_depth_to_surface_1 - sample_depth_upper
        )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("k_2", "hpm_k_2_rep")),
      names_to = "variable",
      values_to = "y"
    ) %>%
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "k_2" ~ "No",
          TRUE ~ "Yes"
        ),
      y_mean = mean(y),
      y_lower = quantile(y, probs = 0.025)[1, ],
      y_upper = quantile(y, probs = 0.975)[1, ],
      x_mean = mean(x),
      x_lower = quantile(x, probs = 0.025)[1, ],
      x_upper = quantile(x, probs = 0.975)[1, ]
    ) %>%
    ggplot(aes(y = y_mean, x = x_mean)) +
    geom_errorbarh(aes(xmin = x_lower, xmax = x_upper), height = 0, color = "grey") +
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), width = 0, color = "grey") +
    geom_smooth(aes(color = variable, fill = variable), formula = y ~ x, method = "lm", se = FALSE) +
    geom_point(aes(fill = variable), shape = 21) +
    facet_grid(taxon_rank_value_pretty ~ id_citation, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA), breaks = equal_breaks(n = 4, s = 0.05, digits = 2)) +
    scale_fill_manual(values = c("black", "grey")) +
    scale_color_manual(values = c("black", "grey")) +
    labs(
      y = expression("Decomposition rate (yr"^{-1}*")"),
      x = "Water table depth below litterbag (cm)"
    ) +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      strip.text.x = ggtext::element_markdown(size = 12),
      strip.text.y = ggtext::element_markdown(angle = 0, hjust = 0, size = 11),
      panel.spacing.y = unit(0.8, "lines")
    ) +
    guides(
      color = guide_legend(title = "Predicted with HPM?", override.aes = list(size = 3)),
      fill = guide_legend(title = "Predicted with HPM?")
    )

  ggsave(
    file_plot,
    plot = res_plot,
    width = 12, height = 7, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}


#' k$_0$ predicted by HPM modifications minus k$_0$ predicted by standard HPM versus water table depth
#'
#' @param x_stan_draws
#'
#' @param hpmd_stan_draws_1
#'
#' @export
hpmd_get_plot_9 <- function(x_stan_draws, hpmd_stan_draws_1) {

  # function to add custom breaks. Modified from https://stackoverflow.com/a/55980394
  equal_breaks <- function(n = 3, s = 0.05, ..., digits){
    function(x){
      # rescaling
      d <- s * diff(range(x)) / (1+2*s)
      round(seq(min(x)+d, max(x)-d, length=n), digits = digits)
    }
  }

  # combine data
  res <-
    x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars) %>%
    dplyr::left_join(
      hpmd_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model =
        model_name |>
        hpmd_model_id_to_name() |>
        factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
    ) |>
    dplyr::left_join(
      hpmd_stan_draws_1 |>
        dplyr::select(dplyr::all_of(c("hpm_k_2_rep", "id_sample"))) |>
        dplyr::rename(
          hpm_k_2_rep_standard = "hpm_k_2_rep"
        ),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      y = hpm_k_2_rep - hpm_k_2_rep_standard
    )

  res_plot <-
    res %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::filter(index_hpm) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start))) %>%
    dplyr::filter(
      taxon_rank_value %in% c("Sphagnum", "Sphagnum angustifolium", "Sphagnum balticum", "Sphagnum cuspidatum", "Sphagnum fuscum", "Sphagnum magellanicum")
    ) %>%
    dplyr::mutate(
      x = layer_water_table_depth_to_surface_1 - sample_depth_upper
    ) %>%
    dplyr::mutate(
      y_mean = mean(y),
      y_lower = quantile(y, probs = 0.025)[1, ],
      y_upper = quantile(y, probs = 0.975)[1, ],
      x_mean = mean(x),
      x_lower = quantile(x, probs = 0.025)[1, ],
      x_upper = quantile(x, probs = 0.975)[1, ]
    ) %>%
    ggplot(aes(y = y_mean, x = x_mean)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_errorbarh(aes(xmin = x_lower, xmax = x_upper), height = 0, color = "grey") +
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), width = 0, color = "grey") +
    # geom_smooth(aes(color = variable, fill = variable), formula = y ~ x, method = "lm", se = FALSE) +
    geom_point(shape = 21) +
    facet_grid(taxon_rank_value_pretty ~ id_model, scales = "free_y") +
    scale_y_continuous(breaks = equal_breaks(n = 4, s = 0.05, digits = 2)) +
    # scale_fill_manual(values = c("black", "grey")) +
    # scale_color_manual(values = c("black", "grey")) +
    labs(
      y = "<i>k<sub>0</sub></i>(HPM modification) - <i>k<sub>0</sub></i>(HPMf) (yr<sup>-1</sup>)",
      x = "Water table depth below litterbag (cm)"
    ) +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      axis.title.y = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      strip.text.y = ggtext::element_markdown(angle = 0, hjust = 0, size = 11),
      panel.spacing.y = unit(0.8, "lines")
    )

  file_plot <- "figures/hpmd_plot_9.pdf"

  ggsave(
    file_plot,
    plot = res_plot,
    width = 6.5, height = 7, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}



#### Tables ####

#' Prior justification table
#'
#' @param x_stan_data
#'
#' @param nsim Positive integer value. The number of simulations to run to
#' compute properties of the prior distributions.
#'
#' @export
hpmd_get_table_prior_justification <- function(x_stan_data, nsim) {

  rvar_inv_logit <- posterior::rfun(binomial()$linkinv)

  # get table part from leaching project
  # res_leaching <- leaching_get_table_prior_justification(x_stan_data = x_stan_data, nsim = nsim)

  group_level_variance_parameters <- c("k_2_p1_p2", "k_2_p2_p2", "k_2_p3_p2", "k_2_p4_p2", "phi_2_p2_p1_p2", "phi_2_p2_p2_p2", "phi_2_p2_p3_p2", "phi_2_p2_p4_p2", "l_2_p1_p2", "l_2_p2_p2", "l_2_p3_p2", "l_2_p4_p2")
  parameters_with_hyperparameters <- c("l_2_p1", "l_2_p2", "l_2_p3", "l_2_p4", "k_2_p1", "k_2_p2", "k_2_p3", "k_2_p4", "phi_2_p2_p1", "phi_2_p2_p2", "phi_2_p2_p3", "phi_2_p2_p4")
  hpm_parameters <- c("layer_total_porosity_1", "layer_minimum_degree_of_saturation_at_surface_1", "layer_water_table_depth_to_surface_1", "hpm_k_2_p1", "m69_p1", "m69_p2", "m68_p1", "m68_p2", "m68_p3_2_p1", "hpm_l_2_p1", "hpm_l_2_p3", "hpm_l_2_p4")

  tibble::tibble(
    parameter_name = c(parameters_with_hyperparameters, "alpha_2_p1", "alpha_2_p2", "alpha_2_p3", "alpha_2_p4", group_level_variance_parameters, hpm_parameters),
    parameter_distribution = c(rep("normal", 16), rep("half-normal", 12), rep("beta", 2), "normal", "gamma", "beta", rep("gamma", 3), rep("normal", 3), "gamma"),
    has_hyperparameter = parameter_name %in% parameters_with_hyperparameters,
    is_hierarchical = parameter_name %in% group_level_variance_parameters,

    parameter_p1 =
      purrr::map_dbl(seq_along(parameter_name), function(i) {
        if(parameter_distribution[[i]] == "half-normal") {
          0
        } else if(parameter_name[[i]] == "layer_water_table_depth_to_surface_1") {
          NA_real_
        } else {
          x_stan_data[[paste0(parameter_name[[i]], "_p1")]] %>%
            unlist() %>%
            unique()
        }
      }),
    parameter_p2 =
      purrr::map_dbl(seq_along(parameter_name), function(i) {
        if(parameter_distribution[[i]] == "half-normal") {
          x_stan_data[[paste0(parameter_name[[i]], "_p1")]] %>%
            unlist() %>%
            unique()
        } else if(has_hyperparameter[[i]]) {
          1
        } else {
          x_stan_data[[paste0(parameter_name[[i]], "_p2")]] %>%
            unlist() %>%
            unique()
        }
      }),
    unit =
      dplyr::case_when(
        stringr::str_detect(parameter_name, "^l_2") ~ "(g g$_\\text{initial}$) (logit scale)",
        stringr::str_detect(parameter_name, "^k_2") ~ "(yr$^{-1}$) (log scale)",
        stringr::str_detect(parameter_name, "^alpha_2") ~ "(-) (log scale)",
        stringr::str_detect(parameter_name, "^phi_2") ~ "(-) (log scale)",
        parameter_name == "layer_total_porosity_1" ~ "L$_\\text{pores}$ L$_\\text{sample}^{-1}$",
        parameter_name == "layer_minimum_degree_of_saturation_at_surface_1" ~
          "L$_\\text{water}$ L$_\\text{pores}^{-1}$",
        parameter_name == "layer_water_table_depth_to_surface_1" ~ "cm",
        parameter_name %in% c("hpm_k_2_p1", "hpm_l_2_p4") ~ "(-)",
        parameter_name == "m69_p1" ~ "L$_\\text{water}$ L$_\\text{pores}^{-1}$",
        parameter_name == "m69_p2" ~ "(-)",
        parameter_name == "m68_p1" ~ "(yr$^{-1}$)",
        parameter_name == "m68_p2" ~ "(cm)",
        parameter_name == "m68_p3_2_p1" ~ "(yr$^{-1}$) (log scale)",
        parameter_name == "hpm_l_2_p1" ~ "(g g$_\\text{initial}^{-1}$) (logit scale)",
        parameter_name == "hpm_l_2_p3" ~ "(g g$_\\text{initial}$ L$_\\text{water}^{-1}$ L$_\\text{pores}$) (logit scale)",
      ),
    hpm_parameter_name =
      dplyr::case_when(
        parameter_name == "m69_p1" ~ "$W_{opt}$",
        parameter_name == "m69_p2" ~ "$c_{1}$",
        parameter_name == "m68_p1" ~ "$f_{min}$",
        parameter_name == "m68_p2" ~ "$c_{2}$",
        parameter_name == "m68_p3_2_p2" ~ "$k_{0,i}$ (log scale)",
        TRUE ~ ""
      ),
    ci95_normal =
      purrr::map(seq_along(parameter_name), function(i) {
        cur_sd <-
          if(has_hyperparameter[[i]]) {
            abs(rnorm(nsim, 0, x_stan_data[[paste0(parameter_name[[i]], "_p2_p1")]]))
          } else {
            parameter_p2[[i]]
          }
        rnorm(
          n = nsim,
          mean = parameter_p1[[i]],
          sd = cur_sd
        ) %>%
          quantile(probs = c(0.025, 0.975), na.rm = TRUE)
      }),
    ci95_beta =
      purrr::map(seq_along(parameter_name), function(i) {
        rbeta(
          n = nsim,
          shape1 = parameter_p1[[i]],
          shape2 = parameter_p2[[i]]
        ) %>%
          quantile(probs = c(0.025, 0.975), na.rm = TRUE)
      }),
    justification =
      dplyr::case_when(
        parameter_name == "l_2_p1" ~ purrr::map_chr(ci95_normal, function(.x) {
          paste0("Assumes an average initial leaching loss across all available litterbag data within (95\\% confidence interval) (", .x %>% binomial(link = "logit")$linkinv() %>% round(3) %>% knitr::combine_words(and = ", "), ") g g$_\\text{initial}^{-1}$")
        }),
        parameter_name == "k_2_p1" ~ purrr::map_chr(ci95_normal, function(.x) {
          paste0("Assumes an average initial decomposition rate across all available litterbag data within (95\\% confidence interval) (", .x %>% exp() %>% round(3) %>% knitr::combine_words(and = ", "), ") yr$^{-1}$")
        }),
        parameter_name == "alpha_2_p1" ~ purrr::map_chr(ci95_normal, function(.x) {
          paste0("Assumes an average $\\alpha$ across all available litterbag data within (95\\% confidence interval) (", .x %>% exp() %>% magrittr::add(1) %>% round(3) %>% knitr::combine_words(and = ", "), ")")
        }),
        hpm_parameter_name != "" & parameter_name != "m68_p3_2_p1" ~ "Centered at the standard value used in the HPM.",
        parameter_name == "layer_total_porosity_1" ~ "(ref:sup-out-d-sdm-all-priors-1-layer-total-porosity)",
        parameter_name == "layer_minimum_degree_of_saturation_at_surface_1" ~ purrr::map_chr(ci95_beta, function(.x) {
          paste0("This parameter comes from the modified Granberg model. The prior distribution assumes a minimum degree of saturation at the surface in different litterbag experiments of (95\\% confidence interval) (", .x %>% round(3) %>% knitr::combine_words(and = ", "), ")")
        }),
        parameter_name == "layer_water_table_depth_to_surface_1" ~ "The average was set to the average water table depths reported in the litterbag studies.",
        parameter_name == "m68_p3_2_p1" ~ purrr::map_chr(ci95_normal, function(.x) {
          paste0("Assumes a maximum potential initial decomposition rate across all species within (95\\% confidence interval) (", .x %>% exp() %>% round(3) %>% knitr::combine_words(and = ", "), ") yr$^{-1}$")
        }),
        parameter_name == "hpm_l_2_p1" ~ purrr::map_chr(ci95_normal, function(.x) {
          paste0("Assumes a maximum possible initial leaching loss across all available litterbag data within (95\\% confidence interval) (", .x %>% rvar_inv_logit() %>% mean() %>% round(3) %>% knitr::combine_words(and = ", "), ") g g$_\\text{initial}^{-1}$")
        }),
        TRUE ~ ""
      )
  ) %>%
    dplyr::mutate(
      parameter_p1 =
        dplyr::case_when(
          parameter_name == "layer_water_table_depth_to_surface_1" ~ "average reported WTD",
          TRUE ~ round(parameter_p1, digits = 2) %>% as.character()
        ),
      parameter_distribution2 =
        paste0(
          parameter_distribution, "(",
          parameter_p1, ", ",
          ifelse(has_hyperparameter, paste0(parameter_name, "_p2"), round(parameter_p2, digits = 2)), ")"
        ) %>%
        stringr::str_replace_all(pattern = "_", replacement = "\\\\_"),
      parameter_name =
        parameter_name %>%
        stringr::str_replace_all(pattern = "_", replacement = "\\\\_")
    ) %>%
    dplyr::select(parameter_name, hpm_parameter_name, unit, parameter_distribution2, justification)

}


