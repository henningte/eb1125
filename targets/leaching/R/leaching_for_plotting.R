#### Helper functions ####

#' Makes pretty taxon_rank_value for plotting
#'
#' @param x
#'
#' @export
leaching_add_taxon_rank_value_pretty <- function(x) {

  res <-
    x %>%
    dplyr::mutate(
      # pretty taxon_rank_value for plotting
      taxon_rank_value_pretty =
        dplyr::case_when(
          taxon_rank_value == "Sphagnum russowii or capillifolium" ~ "S. russowii* and *capillifolium",
          taxon_rank_value == "Sphagnum magellanicum" ~ "S. magellanicum aggr.",
          TRUE ~ taxon_rank_value
        ),
      taxon_rank_value_pretty =
        taxon_rank_value_pretty %>%
        stringr::str_replace("Sphagnum ", replacement = "S. ") %>%
        paste0("*", ., "*") %>%
        stringr::str_replace("^\\*Sphagnum\\*$", replacement = "*Sphagnum* spec.")
    )

  # convert to factor to make sure the ordering is the same

  res <-
    res %>%
    dplyr::mutate(
      taxon_rank_value_pretty =
        taxon_rank_value_pretty %>%
        factor(
          levels =
            res %>%
            dplyr::filter(! duplicated(taxon_rank_value)) %>%
            dplyr::select(taxon_rank_value, taxon_rank_value_pretty) %>%
            dplyr::arrange(taxon_rank_value) %>%
            dplyr::pull(taxon_rank_value_pretty)
        )
    )

  res

}


#' Identifies outliers in predicted average remaining masses
#'
#' @param x_stan_draws A data frame with MCMC draws from a model.
#'
#' @param alpha A numeric value. The significance level to use for the test.
#'
#' @export
leaching_identify_outliers_mass_relative_mass_mu <- function(x_stan_draws, alpha = 0.95) {

  probs = c((1 - alpha) / 2, 1 - (1 - alpha) / 2)
  res <-
    x_stan_draws %>%
    dplyr::mutate(
      .mu_lower = posterior::quantile2(mu, probs = probs[[1]]),
      .mu_upper = posterior::quantile2(mu, probs = probs[[2]]),
      is_outlier = mass_relative_mass > .mu_upper | mass_relative_mass < .mu_lower
    ) %>%
    dplyr::group_by(id_sample_incubation_start) %>%
    dplyr::summarise(
      is_outlier = any(is_outlier),
      .groups = "drop"
    ) %>%
    dplyr::filter(! is.na(is_outlier))

  res

}



#### Plots ####

#' Estimated $k_0$, $l_0$, $\alpha$ for each species and study
#'
#' @param x_stan_draws Data frame with the MCMC draws for the current fit.
#'
#' @param id_fit .
#'
#' @param ... .
#'
#' @export
leaching_plot_estimates_for_species_and_studies <- function(x_stan_draws, id_fit, ...) {

  .dots <- list(...)

  x_stan_draws <-
    x_stan_draws %>%
    dplyr::filter(! is.na(k_2)) %>%
    leaching_add_taxon_rank_value_pretty()

  res <- list()

  # l_2
  if(.dots$has_leaching) {

    res$p1 <-
      x_stan_draws %>%
      dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
      dplyr::mutate(
        id_citation =
          id_citation %>%
          factor() %>%
          as.integer() %>%
          as.factor(),
        l_2 = l_2 * 100,
        l_2_mean = mean(l_2),
        l_2_lower = quantile(l_2, probs = 0.025)[1, ],
        l_2_upper = quantile(l_2, probs = 0.975)[1, ]
      ) %>%
      ggplot(aes(y = l_2_mean, x = id_citation)) +
      # geom_hline(yintercept = c(10, 20, 30, 40, 50, 60), color = "grey", linetype = 2) +
      geom_errorbar(
        aes(ymin = l_2_lower, ymax = l_2_upper),
        width = 0,
        color = "grey",
        position = position_dodge(width = 0.5)
      ) +
      geom_point(
        position = position_dodge(width = 0.5)
      ) +
      facet_grid(~ taxon_rank_value_pretty, switch = "x", scales = "free_x", space = "free") +
      labs(x = "Species", y = "Initial leaching loss (mass-%)") +
      theme(
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y = ggtext::element_markdown(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey"),
        panel.grid.minor.y = element_line(colour = "lightgrey")
      )

  }

  # alpha_2
  if(.dots$has_alpha) {

    res$p2 <-
      x_stan_draws %>%
      dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
      dplyr::mutate(
        id_citation =
          id_citation %>%
          factor() %>%
          as.integer() %>%
          as.factor(),
        l_2 = alpha_2,
        l_2_mean = mean(l_2),
        l_2_lower = quantile(l_2, probs = 0.025)[1, ],
        l_2_upper = quantile(l_2, probs = 0.975)[1, ]
      ) %>%
      ggplot(aes(y = l_2_mean, x = id_citation)) +
      # geom_hline(yintercept = c(1), color = "grey", linetype = 2) +
      geom_errorbar(
        aes(ymin = l_2_lower, ymax = l_2_upper),
        width = 0,
        color = "grey",
        position = position_dodge(width = 0.5)
      ) +
      geom_point(
        position = position_dodge(width = 0.5)
      ) +
      facet_grid(~ taxon_rank_value_pretty, switch = "x", scales = "free_x", space = "free") +
      labs(x = "Species", y = expression(alpha~"(-)")) +
      theme(
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y = ggtext::element_markdown(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey"),
        panel.grid.minor.y = element_line(colour = "lightgrey")
      )

  }

  # k_2
  res$p3 <-
    x_stan_draws %>%
    dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
    dplyr::mutate(
      id_citation =
        id_citation %>%
        factor() %>%
        as.integer() %>%
        as.factor(),
      k_2_mean = mean(k_2),
      k_2_lower = quantile(k_2, probs = 0.025)[1, ],
      k_2_upper = quantile(k_2, probs = 0.975)[1, ]
    ) %>%
    ggplot(aes(y = k_2_mean, x = id_citation)) +
    # geom_hline(yintercept = c(seq(0.1, 0.9, by = 0.1)), color = "grey", linetype = 2) +
    geom_errorbar(
      aes(ymin = k_2_lower, ymax = k_2_upper),
      width = 0,
      color = "grey",
      position = position_dodge(width = 0.5)
    ) +
    geom_point(
      position = position_dodge(width = 0.5)
    ) +
    facet_grid(cols = vars(taxon_rank_value_pretty), switch = "x", scales = "free_x", space = "free") +
    labs(x = "Species", y = expression("Decomposition rate (yr"^{-1}*")")) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown(angle = 90, hjust = 0),
      strip.text.y = ggtext::element_markdown(),
      panel.grid.major.y = element_line(colour = "lightgrey"),
      panel.grid.minor.y = element_line(colour = "lightgrey")
    )


  ## combine
  res_plot <-
    res %>%
    patchwork::wrap_plots(nrow = 1 + .dots$has_leaching + .dots$has_alpha, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    )

  ## metadata for caption
  res_caption <-
    x_stan_draws %>%
    dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
    dplyr::mutate(
      id_citation_number =
        id_citation %>%
        factor() %>%
        as.integer() %>%
        as.factor()
    ) %>%
    dplyr::filter(! duplicated(id_citation_number)) %>%
    dplyr::select(id_citation_number, id_citation) %>%
    dplyr::arrange(as.integer(id_citation_number))

  res_caption <-
    paste0(
      if(.dots$has_leaching & .dots$has_alpha) {
        "Estimated initial leaching losses (a), the parameter controlling a decrease of decomposition rates over time ($\\alpha$) (b), and decomposition rates (c) grouped by species and study for model "
      } else if (.dots$has_leaching) {
        "Estimated initial leaching losses (a), and decomposition rates (c) grouped by species and study for model "
      } else if (.dots$has_alpha) {
        "Estimated parameter controlling a decrease of decomposition rates over time ($\\alpha$) (a), and decomposition rates (b) grouped by species and study for model "
      } else {
        "Estimated  decomposition rates grouped by species and study for model "
      },
      leaching_d_models$model_name[leaching_d_models$id_fit == id_fit], ". ",
      "Points represent averages and error bars 95% confidence intervals. The study is indicated by numbers on the x axis: ",
      res_caption %>%
        dplyr::mutate(
          id_citation_number = paste0("(", id_citation_number, ") "),
          id_citation = paste0(id_citation_number, "@", id_citation)
          ) %>%
        dplyr::pull(id_citation) %>%
        paste(collapse = ", "),
      ". *Sphagnum* spec. are samples that have been identified only to the genus level."
    )

  ## export

  # caption and plot
  file_plot <- paste0("figures/leaching_plot_1_", id_fit, ".pdf")
  file_caption <- paste0("misc_out/leaching_plot_1_", id_fit, "_caption.rds")

  saveRDS(res_caption, file_caption)
  ggsave(
    file_plot,
    plot = res_plot,
    width = 9.5, height = 5 + 2.5 * (.dots$has_leaching + .dots$has_alpha), dpi = 300,
    device = cairo_pdf
  )

  list(
    file_plot = file_plot,
    file_caption = file_caption
  )

}




#' $k_0$ versus $l_0$, each time using the model considering initial leaching losses and the same model not considering initial leaching losses
#'
#' @export
leaching_plot_k0_l0 <- function(x_stan_draws_no_leaching, x_stan_draws_with_leaching, n_experiments_per_species = 5, facet_nrow = 3L, width = 8, height = 9) {

  no_leaching_id_fit <- x_stan_draws_no_leaching$id_fit[[1]]
  with_leaching_id_fit <- x_stan_draws_with_leaching$id_fit[[1]]

  ## combine data on predictions and parameters for plotting
  x_stan_draws_combined <-
    dplyr::bind_rows(
      x_stan_draws_no_leaching %>%
        dplyr::mutate(
          considers_leaching = "No"
        ),
      x_stan_draws_with_leaching %>%
        dplyr::mutate(
          considers_leaching = "Yes"
        )
    ) %>%
    dplyr::filter(! duplicated(paste0(considers_leaching, "_", id_sample_incubation_start))) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::select(-l_2) %>%
    dplyr::left_join(
      x_stan_draws_with_leaching %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
        dplyr::select(id_sample_incubation_start, l_2),
      by = "id_sample_incubation_start"
    )

  # add filter criterion and apply it
  x_stan_draws_combined <-
    dplyr::left_join(
      x_stan_draws_combined,
      x_stan_draws_combined %>%
        dplyr::group_by(considers_leaching, taxon_rank_value_pretty) %>%
        dplyr::summarise(
          index_n_experiments_per_species = length(l_2),
          .groups = "drop"
        ),
      by = c("considers_leaching", "taxon_rank_value_pretty")
    ) %>%
    dplyr::filter(index_n_experiments_per_species >= n_experiments_per_species)

  # add line breaks to long taxon_rank_value_pretty
  levels(x_stan_draws_combined$taxon_rank_value_pretty)[levels(x_stan_draws_combined$taxon_rank_value_pretty) == "*S. russowii* and *capillifolium*"] <- "*S. russowii* and<br>*capillifolium*"
  levels(x_stan_draws_combined$taxon_rank_value_pretty)[levels(x_stan_draws_combined$taxon_rank_value_pretty) == "*S. magellanicum aggr.*"] <- "*S. magellanicum<br>aggr.*"

  ## decomposition rates versus initial leaching losses (all species)

  # plot
  p1 <-
    x_stan_draws_combined %>%
    dplyr::select(id_sample_incubation_start, considers_leaching, k_2, l_2, taxon_rank_value_pretty) %>%
    dplyr::mutate(
      k_2_mean = mean(k_2),
      k_2_lower = quantile(k_2, probs = c(0.025))[1, ],
      k_2_upper = quantile(k_2, probs = c(0.975))[1, ],
      l_2 = mean(l_2 * 100)
    ) %>%
    ggplot(aes(y = k_2_mean, x = l_2)) +
    geom_errorbar(aes(ymin = k_2_lower, ymax = k_2_upper), width = 0, color = "grey") +
    geom_point(aes(fill = considers_leaching), shape = 21) +
    guides(
      fill = guide_legend(title = "Model considers initial leaching losses?", override.aes = list(size = 3))
    ) +
    scale_fill_manual(values = c("black", "grey")) +
    labs(
      y = expression("Decomposition rate (yr"^{-1}*")"),
      x = "Average initial leaching loss (mass-%)"
    ) +
    facet_wrap(~ taxon_rank_value_pretty, nrow = facet_nrow) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown()
    )

  ## standard deviation of decomposition rates versus initial leaching losses (all species)

  # plot
  p2 <-
    x_stan_draws_combined %>%
    dplyr::select(id_sample_incubation_start, considers_leaching, k_2, l_2, taxon_rank_value_pretty) %>%
    dplyr::mutate(
      k_2_sd = sd(k_2),
      l_2 = mean(l_2 * 100)
    ) %>%
    ggplot(aes(y = k_2_sd, x = l_2)) +
    geom_point(aes(fill = considers_leaching), shape = 21) +
    guides(
      fill = guide_legend(title = "Model considers initial leaching losses?", override.aes = list(size = 3))
    ) +
    scale_fill_manual(values = c("black", "grey")) +
    labs(
      y = expression("&sigma;(decomposition rate) (yr<sup>-1</sup>)"),
      x = "Average initial leaching loss (mass-%)"
    ) +
    facet_wrap(~ taxon_rank_value_pretty, nrow = facet_nrow) +
    theme(
      axis.title.y = ggtext::element_markdown(),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = ggtext::element_markdown()
    )

  ## combine
  res_plot <-
    list(p1, p2) %>%
    patchwork::wrap_plots(ncol = 1L, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    ) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ## export
  file_plot <- paste0("figures/leaching_plot_2_fit_", no_leaching_id_fit, "_and_fit_", with_leaching_id_fit, ".pdf")

  ggsave(
    file_plot,
    plot = res_plot,
    width = width, height = height, dpi = 300,
    device = cairo_pdf
  )

  file_plot

}



#' Plots for uncertainty analysis 1
#'
#' @param x_ua_1 Data frame with results of the uncertainty analysis.
#'
#' @param x_stan_draws Data frame. The corresponding MCMC draws. This is
#' needed to get estimates for `l_2`.
#'
#' @export
leaching_get_plot_3 <- function(x_ua_1, x_stan_draws) {

  # add required columns
  x_ua_1 <-
    dplyr::left_join(
      x_ua_1,
      x_stan_draws %>%
        dplyr::select(id_sample, , k_2, l_2, m_rep),
      by = "id_sample"
    )

  ### Histogram of s_k_2
  sdm_ua1_p1 <-
    x_ua_1 %>%
    dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
    ggplot(aes(x = s_k_2_l_2)) +
    geom_histogram(bins = 30) +
    labs(
      y = "Count",
      x = expression(S["j,"*l[0]](k[0]))
    )


  ### Plot: s_k_2 versus duration until first sampling, grouped by leaching losses
  sdm_ua1_p2_data <-
    x_ua_1 %>%
    dplyr::arrange(id_sample_incubation_start, incubation_duration) %>%
    dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
    # dplyr::filter(! id_sample_incubation_start %in% c(24:28, 37:39, 65, 268, 287)) %>%
    dplyr::select(id_sample_incubation_start, s_k_2_l_2, incubation_duration, k_2, l_2, taxon_rank_value) %>%
    dplyr::left_join(
      x_ua_1 %>%
        dplyr::group_by(id_sample_incubation_start) %>%
        dplyr::summarise(
          n = length(l_2),
          m_rep_sd = mean(sd(m_rep)),
          .groups = "drop"
        ),
      by = "id_sample_incubation_start"
    ) %>%
    dplyr::mutate(
      l_2_classes = ggplot2::cut_number(round(mean(l_2 * 100), 0), n = 3),
      k_2_classes = ggplot2::cut_number(round(mean(k_2), 2), n = 3)
    )

  # add correlation coefficients and interval
  sdm_ua1_p2_data <-
    dplyr::left_join(
      sdm_ua1_p2_data,
      sdm_ua1_p2_data %>%
        dplyr::group_by(l_2_classes) %>%
        dplyr::summarise(
          slope_l_2 = {
            res <- summary(lm(s_k_2_l_2 ~ I(incubation_duration/365)))
            paste0("slope = ", round(res$coefficients[2, 1], 2), " (", round(res$coefficients[2, 1] - 2 * res$coefficients[2, 2], 2), ",", round(res$coefficients[2, 1] + 2 * res$coefficients[2, 2], 2), ")")
          },
          .groups = "drop"
        ),
      by = "l_2_classes"
    ) %>%
    dplyr::left_join(
      sdm_ua1_p2_data %>%
        dplyr::group_by(k_2_classes) %>%
        dplyr::summarise(
          slope_k_2 = {
            res <- summary(lm(s_k_2_l_2 ~ I(incubation_duration/365)))
            paste0("slope = ", round(res$coefficients[2, 1], 2), " (", round(res$coefficients[2, 1] - 2 * res$coefficients[2, 2], 2), ",", round(res$coefficients[2, 1] + 2 * res$coefficients[2, 2], 2), ")")
          },
          .groups = "drop"
        ),
      by = "k_2_classes"
    )

  sdm_ua1_p2 <-
    sdm_ua1_p2_data %>%
    ggplot(aes(y = s_k_2_l_2, x = incubation_duration)) +
    # geom_point(aes(size = n, color = m_rep_sd)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey50") +
    geom_text(
      aes(label = slope_l_2, y = 0.9),
      x = 20, hjust = 0
    ) +
    labs(
      y = expression(S["j,"*l[0]](k[0])),
      x = "Incubation duration (d)"
    ) +
    coord_cartesian(xlim = c(0, NA)) +
    facet_wrap(~ l_2_classes, nrow = 1) +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12)
    )



  ### Plot: s_k_2 versus duration until first sampling, grouped by decomposition rates
  sdm_ua1_p3 <-
    sdm_ua1_p2_data %>%
    ggplot(aes(y = s_k_2_l_2, x = incubation_duration)) +
    # geom_point(aes(size = n, color = m_rep_sd)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey50") +
    geom_text(
      aes(label = slope_k_2, y = 0.9),
      x = 20, hjust = 0
    ) +
    labs(
      y = expression(S["j,"*l[0]](k[0])),
      x = "Incubation duration (d)"
    ) +
    coord_cartesian(xlim = c(0, NA)) +
    facet_wrap(~ k_2_classes, nrow = 1) +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12)
    )


  ### combine
  sdm_ua1_p1_p2_p3 <-
    list(
      list(patchwork::plot_spacer(), sdm_ua1_p1, patchwork::plot_spacer()) %>%
        patchwork::wrap_plots(nrow = 1L, byrow = TRUE, widths = c(1, 2, 1)),
      list(sdm_ua1_p2, sdm_ua1_p3) %>%
        patchwork::wrap_plots(nrow = 2L, byrow = TRUE)
    ) %>%
    patchwork::wrap_plots(nrow = 2L, byrow = TRUE, heights = c(1, 3))+
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    )

  ## export
  res_file <- "figures/leaching_plot_3.pdf"
  ggsave(
    res_file,
    plot = sdm_ua1_p1_p2_p3,
    width = 7, height = 8, dpi = 300,
    device = cairo_pdf
  )

  res_file

}



#' Plot the difference of parameter estimates between two models
#'
#' @param x_stan_draws
#'
#' @param y_stan_draws
#'
#' @export
leaching_get_plot_4_individual <- function(x_stan_draws, y_stan_draws) {

  x_id_fit <- x_stan_draws$id_fit[[1]]
  y_id_fit <- y_stan_draws$id_fit[[1]]
  x_model_name <- leaching_d_models$model_name[leaching_d_models$id_fit == x_id_fit]
  y_model_name <- leaching_d_models$model_name[leaching_d_models$id_fit == y_id_fit]
  res_file <- paste0("figures/leaching_plot_4_fit_", x_id_fit, "_vs_fit_", y_id_fit, ".pdf")

  # data
  res <-
    dplyr::left_join(
      x_stan_draws %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)),
      y_stan_draws %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
        dplyr::select(l_2, k_2, id_sample) %>%
        dplyr::rename(
          l_2_fit_y = "l_2",
          k_2_fit_y = "k_2"
        ),
      by = "id_sample"
    ) %>%
    dplyr::mutate(
      l_2_diff = (l_2 - l_2_fit_y) * 100,
      k_2_diff = k_2 - k_2_fit_y
    ) %>%
    dplyr::filter(! is.na(k_2_diff))

  # initial leaching
  p1 <-
    res %>%
    ggplot(aes(ydist = l_2_diff, x = seq_len(nrow(res)))) +
    ggdist::stat_lineribbon() +
    scale_fill_brewer() +
    geom_hline(yintercept = 0, color = "grey50") +
    guides(
      fill = guide_legend(title = "Confidence level")
    ) +
    labs(
      y = paste0("<i>l</i><sub>0, model ", x_model_name, "</sub> - <i>l</i><sub>0, model ", y_model_name, "</sub> (%)"),
      x = "Sample number (-)"
    ) +
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.text = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      legend.position = "bottom"
    )


  # decomposition rate
  p2 <-
    res %>%
    ggplot(aes(ydist = k_2_diff, x = seq_len(nrow(res)))) +
    ggdist::stat_lineribbon() +
    scale_fill_brewer() +
    geom_hline(yintercept = 0, color = "grey50") +
    guides(
      fill = guide_legend(title = "Confidence level")
    ) +
    labs(
      y = paste0("<i>k</i><sub>0, model ", x_model_name, "</sub> - <i>k</i><sub>0, model ", y_model_name, "</sub> (yr<sup>-1</sup>)"),
      x = "Sample number (-)"
    ) +
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.text = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      legend.position = "bottom"
    )


  # combine
  res_combined <-
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
    res_file,
    plot = res_combined,
    width = 9, height = 6, dpi = 300,
    device = cairo_pdf
  )

  res_file

}


#' Plots littebag trajectories of outliers
#'
#' @param x_stan_draws
#'
#' @param alpha
#'
#' @param id_fit
#'
#' @export
leaching_get_plot_outlier_trajectories <- function(x_stan_draws, alpha, width, height) {

  id_fit <- x_stan_draws$id_fit[[1]]

  # identify outliers
  res <-
    dplyr::left_join(
      x_stan_draws,
      x_stan_draws %>%
        leaching_identify_outliers_mass_relative_mass_mu(alpha = alpha),
      by = "id_sample_incubation_start"
    ) %>%
    dplyr::filter(! is.na(m_rep)) %>%
    dplyr::filter(is_outlier) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::mutate(
      id_citation_number =
        id_citation %>%
        factor() %>%
        as.integer() %>%
        as.factor()
    )

  # plot
  res_plot <-
    res %>%
    ggplot(aes(x = incubation_duration/365)) +
    ggdist::stat_ribbon(aes(ydist = m_rep * 100)) +
    geom_path(aes(y = mass_relative_mass * 100), color = "red") +
    scale_fill_brewer() +
    facet_wrap(~ id_sample_incubation_start + taxon_rank_value_pretty + id_citation_number) +
    labs(
      y = "Remaining mass (%)",
      x = "Incubation duration (yr)"
    ) +
    labs(
      fill = guide_legend(title = "Confidence level")
    ) +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      legend.position = "bottom"
    )

  ## metadata for caption
  res_caption <-
    res %>%
    dplyr::filter(! duplicated(id_citation_number)) %>%
    dplyr::select(id_citation_number, id_citation) %>%
    dplyr::arrange(as.integer(id_citation_number))

  ## caption
  res_caption <-
    paste0(
      "Remaining masses predicted by model ", leaching_d_models$model_name[leaching_d_models$id_fit == id_fit],
      " (shaded area) or as reported on average in available litterbag studies (red lines) during the litterbag experiments for outlier litterbag replicates. There is a panel for each litterbag replicate, identified by `id_sample_incubation_start` in the Peatland Decomposition Database (first row), the species (second row), and the study the data are from (third) row. ",
      "The studies are: ", res_caption |> dplyr::mutate(id_citation_number = paste0("(", id_citation_number, ") "), id_citation = paste0(id_citation_number, "@", id_citation)) |> dplyr::pull(id_citation) |> paste(collapse = ", "),
      ". Outlier litterbag replicates are defined as those replicates where the average measured remaining mass significantly different from the average predicted remaining mass ($\\alpha = ", alpha, "$)."
    )

  ## export
  file_plot <- paste0("figures/leaching_plot_outlier_1_", id_fit, ".pdf")
  file_caption <- paste0("misc_out/leaching_plot_outlier_1_", id_fit, "_caption.rds")

  saveRDS(res_caption, file_caption)
  ggsave(
    file_plot,
    plot = res_plot,
    width = width, height = height, dpi = 300,
    device = cairo_pdf
  )

  list(
    plot = file_plot,
    caption = file_caption
  )

}



#' Plot to illustrate that different values of l_0, k_0, and alpha can produce similar mass-time trajectories in littebag experiments
#'
#' @export
leaching_get_simulation_interpretation_alpha_1 <- function() {

  leaching_simulation_interpretation_alpha_1_data <-
    dplyr::bind_rows(
      tibble::tibble(
        correct_interpretation = TRUE,
        t = seq(0, 10, length.out = 100),
        is_extrapolation = t > 5,
        l0 = 0.1,
        k0 = 0.1,
        alpha = 2,
        m = {
          res <- (1 - l0) / (1 + (alpha - 1) * k0 * t)^(1/(alpha - 1))
          res[t <= 0] <- 1
          res
        }
      ),
      tibble::tibble(
        correct_interpretation = FALSE,
        t = seq(0, 10, length.out = 100),
        is_extrapolation = t > 5,
        l0 = 0,
        k0 = 0.5,
        alpha = 6,
        m = {
          (1 - l0) / (1 + (alpha - 1) * k0 * t)^(1/(alpha - 1))
        }
      )
    ) |>
    dplyr::mutate(
      variable = paste0("<i>l<sub>0</sub></i> = ", l0 * 100, "%, <i>k<sub>0</sub></i> = ", k0, " yr<sup>-1</sup>, &alpha; = ", round(alpha, 2))
    )

  res_plot <-
    leaching_simulation_interpretation_alpha_1_data |>
    ggplot(aes(y = m * 100, x = t, color = variable, linetype = is_extrapolation)) +
    geom_path() +
    labs(
      y = "Remaining mass (%)",
      x = "Incubation duration (yr)"
    ) +
    scale_color_manual(values = c("grey", "black")) +
    guides(
      color = guide_legend(title = "", nrow = 2, byrow = TRUE),
      linetype = "none"
    ) +
    theme(
      legend.text = ggtext::element_markdown(),
      legend.position = "bottom"
    )

  file_plot <- paste0("figures/leaching_simulation_interpretation_alpha_1.pdf")

  ggsave(
    file_plot,
    plot = res_plot,
    width = 5.5, height = 3.8, dpi = 300,
    device = cairo_pdf
  )

  list(
    leaching_simulation_interpretation_alpha_1_data = leaching_simulation_interpretation_alpha_1_data,
    file_plot = file_plot
  )

}



#' Plots estimated and true simulated parameter values for `leaching_stan_draws_13`
#'
#' @param x_stan_draws
#'
#' @export
leaching_get_plot_5 <- function(x_stan_draws) {

  p1 <-
    x_stan_draws %>%
    dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
    dplyr::mutate(
      l_2_sim = l_2_sim * 100,
      l_2 = l_2 * 100
    ) %>%
    # dplyr::filter(n_samples_per_sampling_date == 5) |>
    ggplot(aes(ydist = k_2_sim - k_2, x = l_2_sim, color = as.factor(round(alpha_2_sim, 2)))) +
    ggdist::stat_pointinterval(.width = 0.95, interval_size = 0.5) +
    geom_hline(yintercept = 0, color = "grey50") +
    scale_color_manual(values = c("black", "grey")) +
    facet_grid(k_2_sim ~ incubation_duration_design) +
    labs(
      y = "<i>k</i><sub>0, true</sub> - <i>k</i><sub>0, model 1-6</sub> (yr<sup>-1</sup>)",
      x = "<i>l</i><sub>0, true</sub> (mass-%)"
    ) +
    guides(
      color = guide_legend(title = "<i>&alpha;</i><sub>true</sub>")
    ) +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.position = "bottom"
    )

  p2 <-
    x_stan_draws |>
    dplyr::filter(! duplicated(id_sample_incubation_start)) |>
    dplyr::mutate(
      l_2_sim = l_2_sim * 100,
      l_2 = l_2 * 100
    ) %>%
    # dplyr::filter(n_samples_per_sampling_date == 5) |>
    ggplot(aes(ydist = l_2_sim - l_2, x = l_2_sim, color = as.factor(round(alpha_2_sim, 2)))) +
    ggdist::stat_pointinterval(.width = 0.95, interval_size = 0.5) +
    geom_hline(yintercept = 0, color = "grey50") +
    scale_color_manual(values = c("black", "grey")) +
    facet_grid(k_2_sim ~ incubation_duration_design) +
    labs(
      y = "<i>l</i><sub>0, true</sub> - <i>l</i><sub>0, model 1-6</sub> (mass-%)",
      x = "<i>l</i><sub>0, true</sub> (mass-%)"
    ) +
    guides(
      color = guide_legend(title = "<i>&alpha;</i><sub>true</sub>")
    ) +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.position = "bottom"
    )

  p3 <-
    x_stan_draws |>
    dplyr::filter(! duplicated(id_sample_incubation_start)) |>
    dplyr::mutate(
      l_2_sim = l_2_sim * 100,
      l_2 = l_2 * 100
    ) %>%
    # dplyr::filter(n_samples_per_sampling_date == 5) |>
    ggplot(aes(ydist = alpha_2_sim - alpha_2, x = l_2_sim, color = as.factor(round(alpha_2_sim, 2)))) +
    ggdist::stat_pointinterval(.width = 0.95, interval_size = 0.5) +
    geom_hline(yintercept = 0, color = "grey50") +
    scale_color_manual(values = c("black", "grey")) +
    facet_grid(k_2_sim ~ incubation_duration_design) +
    labs(
      y = "<i>&alpha;</i><sub>true</sub> - <i>&alpha;</i><sub>model 1-6</sub> (-)",
      x = "<i>l</i><sub>0, true</sub> (mass-%)"
    ) +
    guides(
      color = guide_legend(title = "<i>&alpha;</i><sub>true</sub>")
    ) +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.position = "bottom"
    )

  file_plot <- c(paste0("figures/leaching_plot_5_plot_", 1:3, ".pdf"))

  for(i in seq_along(file_plot)) {
    ggsave(
      file_plot[[i]],
      plot = get(paste0("p", i)),
      width = 7, height = 4.5, dpi = 300,
      device = cairo_pdf
    )
  }

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
leaching_get_table_prior_justification <- function(x_stan_data, nsim) {

  group_level_variance_parameters <- c("k_2_p1_p2", "k_2_p2_p2", "k_2_p3_p2", "k_2_p4_p2", "phi_2_p2_p1_p2", "phi_2_p2_p2_p2", "phi_2_p2_p3_p2", "phi_2_p2_p4_p2", "l_2_p1_p2", "l_2_p2_p2", "l_2_p3_p2", "l_2_p4_p2")
  parameters_with_hyperparameters <- c("l_2_p1", "l_2_p2", "l_2_p3", "l_2_p4", "k_2_p1", "k_2_p2", "k_2_p3", "k_2_p4", "phi_2_p2_p1", "phi_2_p2_p2", "phi_2_p2_p3", "phi_2_p2_p4")


    tibble::tibble(
      parameter_name = c(parameters_with_hyperparameters, "alpha_2_p1", "alpha_2_p2", "alpha_2_p3", "alpha_2_p4", group_level_variance_parameters),
      parameter_distribution = c(rep("normal", 16), rep("half-normal", 12)),
      has_hyperparameter = parameter_name %in% parameters_with_hyperparameters,
      is_hierarchical = parameter_name %in% group_level_variance_parameters,
      parameter_p1 =
        purrr::map_dbl(seq_along(parameter_name), function(i) {
          if(parameter_distribution[[i]] == "half-normal") {
            0
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
          stringr::str_detect(parameter_name, "^phi_2") ~ "(-) (log scale)"
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
          TRUE ~ ""
        )
    ) %>%
    dplyr::mutate(
      parameter_distribution2 =
        paste0(
          parameter_distribution, "(",
          round(parameter_p1, digits = 2), ", ",
          ifelse(has_hyperparameter, paste0(parameter_name, "_p2"), round(parameter_p2, digits = 2)), ")"
        ) %>%
        stringr::str_replace_all(pattern = "_", replacement = "\\\\_"),
      parameter_name =
        parameter_name %>%
        stringr::str_replace_all(pattern = "_", replacement = "\\\\_")
    ) %>%
    dplyr::select(parameter_name, unit, parameter_distribution2, justification)

}


#' Table with $k_0$, $l_0$, $\alpha$ estimates for all models
#'
#' @param x_stan_draws List of data frames with MCMC draws from the models.
#'
#' @export
leaching_get_table_estimates <- function(x_stan_draws) {

  # combine
  res <-
    x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars) %>%
    dplyr::left_join(
      leaching_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model = paste0("model ", model_name)
    ) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::select(id_sample_incubation_start, taxon_rank_value, taxon_rank_value_pretty, id_citation, id_model, l_2, k_2, alpha_2) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", id_sample_incubation_start)))

  # Adjust columns for printing the tables
  res <-
    res %>%
    dplyr::mutate(
      id_citation =
        paste0("(ref:", stringr::str_remove(id_citation, pattern = "\\."), ")"),
      taxon_rank_value_pretty =
        taxon_rank_value_pretty %>%
        stringr::str_replace_all(pattern = "\\*", replacement = "$") %>%
        stringr::str_replace_all(pattern = "S\\. ", replacement = "S.~") %>%
        stringr::str_replace_all(pattern = "magellanicum aggr", replacement = "magellanicum~aggr")
    ) %>%
    dplyr::select(-taxon_rank_value) %>%
    dplyr::arrange(id_model, taxon_rank_value_pretty)

  list(
    k0 =
      res %>%
      dplyr::select(-l_2, -alpha_2) %>%
      tidyr::pivot_wider(
        names_from = "id_model",
        values_from = c(k_2)
      ) %>%
      dplyr::group_by(taxon_rank_value_pretty, id_citation) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(c("model")),
          function(.x) {
            if(all(is.na(.x))) {
              return("")
            }
            res <- range(mean(.x, na.rm = TRUE))
            paste(round(res, 2), collapse = ", ")
          }
        ),
        n = length(id_sample_incubation_start),
        .groups = "drop"
      ),
    l0 =
      res %>%
      dplyr::select(-k_2, -alpha_2) %>%
      tidyr::pivot_wider(
        names_from = "id_model",
        values_from = c(l_2)
      ) %>%
      dplyr::group_by(taxon_rank_value_pretty, id_citation) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(c("model")),
          function(.x) {
            if(all(is.na(.x))) {
              return("")
            }
            res <- range(mean(.x, na.rm = TRUE))
            if(stringr::str_detect(dplyr::cur_column(), "^model")) {
              res <- res * 100
            }
            paste(round(res, 2), collapse = ", ")
          }
        ),
        n = length(id_sample_incubation_start),
        .groups = "drop"
      ),
    alpha =
      res %>%
      dplyr::select(-l_2, -k_2) %>%
      tidyr::pivot_wider(
        names_from = "id_model",
        values_from = c(alpha_2)
      ) %>%
      dplyr::group_by(taxon_rank_value_pretty, id_citation) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(c("model")),
          function(.x) {
            if(all(is.na(.x))) {
              return("")
            }
            res <- range(mean(.x, na.rm = TRUE))
            paste(round(res, 2), collapse = ", ")
          }
        ),
        n = length(id_sample_incubation_start),
        .groups = "drop"
      )
  ) %>%
    purrr::map(function(.x) {
      .x %>%
        dplyr::relocate(n, .after = "id_citation")
    })

}



#' Table with average estimates for each species
#'
#' @param x_stan_draws List of data frames with extracted MCMC draws.
#'
#' @export
leaching_get_table_estimates_species <- function(x_stan_draws) {

  # combine
  x_stan_draws %>%
    purrr::map_dfr(readRDS_rvars) %>%
    dplyr::left_join(
      leaching_d_models %>%
        dplyr::select(id_fit, model_name),
      by = "id_fit"
    ) %>%
    dplyr::mutate(
      id_model = paste0("model ", model_name)
    ) %>%
    leaching_add_taxon_rank_value_pretty() %>%
    dplyr::select(taxon_rank_value, taxon_rank_value_pretty, id_model, l_2_species, k_2_species, alpha_2_species) %>%
    dplyr::filter(! duplicated(paste0(id_model, "_", taxon_rank_value))) %>%
    dplyr::select(-taxon_rank_value) %>%
    dplyr::mutate(
      taxon_rank_value_pretty =
        taxon_rank_value_pretty %>%
        stringr::str_replace_all(pattern = "\\*", replacement = "$") %>%
        stringr::str_replace_all(pattern = "S\\. ", replacement = "S.~") %>%
        stringr::str_replace_all(pattern = "magellanicum aggr", replacement = "magellanicum~aggr"),
      l_2_species = l_2_species * 100,
      dplyr::across(
        dplyr::all_of(c("l_2_species", "k_2_species", "alpha_2_species")),
        function(.x) {
          cur_digits <-
            switch(
              cur_column(),
              "k_2_species" = 2,
              "l_2_species" = 1,
              "alpha_2_species" = 1
            )

          ifelse(
            is.na(mean(.x)),
            "",
            paste0(round(mean(.x, na.rm = TRUE), digits = cur_digits), " (", round(posterior::quantile2(.x, probs = c(0.025), na.rm = TRUE), digits = cur_digits), ", ", round(posterior::quantile2(.x, probs = c(0.975), na.rm = TRUE), digits = cur_digits), ")")
          )

        }
      )
    ) %>%
    tidyr::pivot_wider(
      names_from = "id_model",
      values_from = dplyr::all_of(c("l_2_species", "k_2_species", "alpha_2_species"))
    )

}
