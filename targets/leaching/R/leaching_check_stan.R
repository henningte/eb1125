#### Plots ####

#' Prior and posterior predictive checks (density plots)
#'
#' @param x Data frame with a numeric column for measured values and a `rvars`
#' column with predicted values
#'
#' @param y_name Character value. Name of the column with measured values.
#'
#' @param yrep_name Character value. Name of the column with predicted values.
#'
#' @param x_lab Character value. X axis title.
#'
#' @param x_scale Numeric value. Used to scale the x axis. Defaults to 1.
#'
#' @param model_name Character value. Name of the model. This is put in the title.
#'
#' @param add_prefix_model Logical vlaue. If `TRUE`, "model " will be added as prefix to the title of the plot.
#'
#' @export
leaching_check_stan_plot_ppc <- function(x, y_name, yrep_name, x_lab, x_scale = 1, model_name, add_prefix_model = TRUE) {

  bayesplot::ppc_dens_overlay(
    y =
      x %>%
      dplyr::filter(! is.na(!!rlang::sym(yrep_name)) & ! is.na(!!rlang::sym(y_name))) %>%
      dplyr::pull(dplyr::all_of(y_name)),
    yrep =
      x %>%
      dplyr::filter(! is.na(!!rlang::sym(yrep_name)) & ! is.na(!!rlang::sym(y_name))) %>%
      dplyr::select(dplyr::all_of(yrep_name)) %>%
      dplyr::pull(dplyr::all_of(yrep_name)) %>%
      posterior::draws_of() %>%
      as.data.frame() %>%
      dplyr::slice_sample(n = 100, replace = FALSE) %>%
      as.matrix()
  ) +
    ggplot2::labs(
      title =
        if(add_prefix_model) {
          paste0("model ", model_name)
        } else {
          paste0(model_name)
        },
      y = "Density (-)",
      x = x_lab
    ) +
    ggplot2::scale_x_continuous(
      labels = function(.x) .x * x_scale
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggtext::element_markdown()
      )

}



#' Measured versus fitted values
#'
#' @param x Data frame with a numeric column for measured values ("y") and a
#' `rvars` column with predicted values ("yrep").
#'
#' @export
leaching_check_stan_plot_y_vs_yrep <- function(x, y_name, yrep_name, x_scale = 1, axis_title_unit, model_name, use_coord_fixed = TRUE) {

  res <-
    x %>%
    dplyr::select(!!rlang::sym(y_name), !!rlang::sym(yrep_name)) %>%
    dplyr::rename(
      y = !!y_name,
      yrep = !!yrep_name
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("y", "yrep")),
        function(.x) .x * x_scale
      )
    ) %>%
    dplyr::select(dplyr::all_of(c("y", "yrep"))) %>%
    as.data.frame() %>%
    ggplot2::ggplot(ggplot2::aes(y = y, xdist = yrep)) +
    ggdist::stat_pointinterval(.width = 0.95, point_interval = ggdist::mean_qi, interval_size = 1, interval_color = "grey") +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey50") +
    ggplot2::labs(
      title = paste0("model ", model_name),
      y = paste0("Measured (", axis_title_unit, ")"),
      x = paste0("Predicted (", axis_title_unit, ")")
    ) +
    ggplot2::theme(
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown()
    )

  if(use_coord_fixed) {
    res <-
      res + ggplot2::coord_fixed()
  }

  res

}




#### Combine individual plots into one plot ####

#' Combining `leaching_check_stan_plot_ppc` for multiple models
#'
#' @param x List of targets from `leaching_stan_draws`.
#'
#' @inheritParams leaching_check_stan_plot_ppc
#'
#' @param file Character value. File where to save the plot (must be pdf).
#'
#' @param width,height Numeric values, the width and height of the pdf.
#'
#' @export
leaching_check_stan_plot_ppc_combined <- function(x, y_name, yrep_name, x_lab, x_scale = 1, file, width, height, add_prefix_model = TRUE) {

  res <-
    x %>%
    purrr::map(function(.x2) {

      leaching_check_stan_plot_ppc(
        x = .x2,
        y_name = y_name,
        yrep_name = yrep_name,
        x_lab = x_lab,
        x_scale = x_scale,
        model_name =
          leaching_d_models %>%
          dplyr::filter(id_fit == .x2$id_fit[[1]]) %>%
          dplyr::pull(model_name),
        add_prefix_model = add_prefix_model
      )

    }) %>%
    patchwork::wrap_plots(ncol = 4L, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    )

  ggplot2::ggsave(
    filename = file,
    plot = res,
    width = width, height = height, dpi = 300,
    device = cairo_pdf
  )

  file

}


#' Combining `leaching_check_stan_plot_y_vs_yrep` for multiple models
#'
#' @param x List of targets from `leaching_stan_draws`.
#'
#' @inheritParams leaching_check_stan_plot_y_vs_yrep
#'
#' @param file Character value. File where to save the plot (must be pdf).
#'
#' @param width,height Numeric values, the width and height of the pdf.
#'
#' @export
leaching_check_stan_plot_y_vs_yrep_combined <- function(x, y_name, yrep_name, x_scale = 1, axis_title_unit, use_coord_fixed = TRUE, file, width, height) {

  res <-
    x %>%
    purrr::map(function(.x2) {

      leaching_check_stan_plot_y_vs_yrep(
        x = .x2,
        y_name = y_name,
        yrep_name = yrep_name,
        x_scale = x_scale,
        axis_title_unit = axis_title_unit,
        use_coord_fixed = use_coord_fixed,
        model_name =
          leaching_d_models %>%
          dplyr::filter(id_fit == .x2$id_fit[[1]]) %>%
          dplyr::pull(model_name)
      )

    }) %>%
    patchwork::wrap_plots(ncol = 4L, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    )

  ggplot2::ggsave(
    filename = file,
    plot = res,
    width = width, height = height, dpi = 300,
    device = cairo_pdf
  )

  file

}


#### Tables ####

#' MCMC diagnostics for one model
#'
#' @export
leaching_check_stan_diagnostics <- function(x_stan_model, id_fit) {

  # get diagnostics. ---note: Use monitor() here because this computes the improved rhat and ess diagnostics and MCSE
  res <-
    rstan::monitor(
      rstan::extract(
        x_stan_model,
        permute = F,
        inc_warmup = F
      ),
      print = FALSE,
      warmup = 0
    )

  res %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      par= rownames(res),
      id_fit = id_fit
    ) %>%
    dplyr::relocate(dplyr::all_of(c("id_fit", "par")), .before = dplyr::everything())

}



#' Summarizes key MCMC diagnostics across parameters for one model
#'
#' @param leaching_table_mcmc_diagnostics Data frame. Output from
#' `leaching_check_stan_diagnostics()`.
#'
#' @export
leaching_check_stan_diagnostics_summary <- function(x_stan_model, leaching_table_mcmc_diagnostics, id_fit) {

  tibble::tibble(
    id_fit = id_fit,
    num_divergent = rstan::get_num_divergent(x_stan_model),
    max_rhat_improved = max(leaching_table_mcmc_diagnostics$Rhat),
    min_ess_bulk = min(leaching_table_mcmc_diagnostics$Bulk_ESS),
    min_ess_tail = min(leaching_table_mcmc_diagnostics$Tail_ESS)
  )

}



#### Combine all ####

#' Wrapper function which combines all checks
#'
#' @export
leaching_stan_get_all_checks <- function(
    prior_models,
    posterior_models,
    x_stan_models
) {

  prior_models <- purrr::map(prior_models, readRDS_rvars)
  posterior_models <- purrr::map(posterior_models, readRDS_rvars)

  res <-
    tibble::lst(
      plots =
        list(
          # m_rep: prior
          leaching_plot_ppc_prior_m =
            prior_models %>%
            leaching_check_stan_plot_ppc_combined(
              y_name = "mass_relative_mass",
              yrep_name = "m_rep",
              x_lab = "Remaining mass (mass-%)",
              x_scale = 100,
              file = "figures/leaching_plot_ppc_prior_m.pdf",
              width = 10,
              height = 3
            ),
          # phi: prior
          leaching_plot_ppc_prior_phi =
            prior_models %>%
            leaching_check_stan_plot_ppc_combined(
              y_name = "mass_relative_mass_precision",
              yrep_name = "phi",
              x_lab = "&Phi; (-)",
              x_scale = 1,
              file = "figures/leaching_plot_ppc_prior_phi.pdf",
              width = 10,
              height = 3
            ),
          # m_rep: posterior
          leaching_plot_ppc_posterior_m =
            posterior_models %>%
            leaching_check_stan_plot_ppc_combined(
              y_name = "mass_relative_mass",
              yrep_name = "m_rep",
              x_lab = "Remaining mass (mass-%)",
              x_scale = 100,
              file = "figures/leaching_plot_ppc_posterior_m.pdf",
              width = 10,
              height = 6
            ),
          # phi: posterior
          leaching_plot_ppc_posterior_phi =
            posterior_models %>%
            leaching_check_stan_plot_ppc_combined(
              y_name = "mass_relative_mass_precision",
              yrep_name = "phi",
              x_lab = "&Phi; (-)",
              x_scale = 1,
              file = "figures/leaching_plot_ppc_posterior_phi.pdf",
              width = 10,
              height = 6
            ),
          # mrep: y vs yrep
          leaching_plot_y_vs_yrep_m =
            posterior_models %>%
            leaching_check_stan_plot_y_vs_yrep_combined(
              y_name = "mass_relative_mass",
              yrep_name = "m_rep",
              x_scale = 100,
              axis_title_unit = "mass-%",
              use_coord_fixed = FALSE,
              file = "figures/leaching_plot_y_vs_yrep_m.pdf",
              width = 11,
              height = 6
            ),
          # phi: y vs yrep
          leaching_plot_y_vs_yrep_phi =
            posterior_models %>%
            leaching_check_stan_plot_y_vs_yrep_combined(
              y_name = "mass_relative_mass_precision",
              yrep_name = "phi",
              x_scale = 1,
              axis_title_unit = "-",
              use_coord_fixed = FALSE,
              file = "figures/leaching_plot_y_vs_yrep_phi.pdf",
              width = 11,
              height = 6
            )
        ),
      tables =
        tibble::lst(
          leaching_table_mcmc_diagnostics =
            purrr::map2(x_stan_models, leaching_d_models$id_fit[! leaching_d_models$prior_only], function(.x, .y) {
              leaching_check_stan_diagnostics(
                x_stan_model = .x,
                id_fit = .y
              )
            }),
          leaching_table_mcmc_diagnostics_summary =
            purrr::map2_dfr(x_stan_models, seq_along(leaching_d_models$id_fit[! leaching_d_models$prior_only]), function(.x, .y) {
              leaching_check_stan_diagnostics_summary(
                x_stan_model = .x,
                leaching_table_mcmc_diagnostics = leaching_table_mcmc_diagnostics[[.y]],
                id_fit = leaching_d_models$id_fit[! leaching_d_models$prior_only][[.y]]
              )
            }),
          leaching_table_mcse_summary =
            purrr::map_dfr(leaching_table_mcmc_diagnostics, function(.x) {

              .x %>%
                dplyr::filter(stringr::str_detect(par, "^l_2\\[") | stringr::str_detect(par, "^k_2\\[") | stringr::str_detect(par, "^alpha_2\\[") | stringr::str_detect(par, "^m_rep\\[")) %>%
                dplyr::mutate(
                  base_par =
                    dplyr::case_when(
                      stringr::str_detect(par, "^l_2\\[") ~ "l_2",
                      stringr::str_detect(par, "^k_2\\[") ~ "k_2",
                      stringr::str_detect(par, "^alpha_2\\[") ~ "alpha_2",
                      stringr::str_detect(par, "^m_rep\\[") ~ "m_rep"
                    )
                ) %>%
                dplyr::group_by(base_par) %>%
                dplyr::summarise(
                  id_fit = unique(id_fit),
                  dplyr::across(dplyr::starts_with("MCSE_"), max, na.rm = TRUE),
                  .groups = "drop"
                )
            })
        )
    )

  res

}



#### Simulations ####

#' Checks whether Bayesian posterior intervals contain the true value
#'
#' @param x_stan_draws
#'
#' @export
leaching_check_coverage_simulation <- function(x_stan_draws, alpha = 0.95) {

  probs <- c((1 - alpha) / 2, 1 - (1 - alpha) / 2)

  x_stan_draws |>
    dplyr::mutate(
      k_2_covered =
        {
          res <- posterior::quantile2(k_2, probs = probs)
          res[1, ] <= k_2_sim & res[2, ] >= k_2_sim
        },
      l_2_covered =
        {
          res <- posterior::quantile2(l_2, probs = probs)
          res[1, ] <= l_2_sim & res[2, ] >= l_2_sim
        },
      alpha_2_covered =
        {
          res <- posterior::quantile2(alpha_2, probs = probs)
          res[1, ] <= alpha_2_sim & res[2, ] >= alpha_2_sim
        }
    )

}



#### priorsense ####

#' Performs a prior and likelihood sensitivity analysis via power scaling as implemented in the priorsense package
#'
#' @param x_stan_fit
#'
#' @export
leaching_do_priorsense_analysis <- function(x_stan_fit, id_fit, width, height) {

  # variables for which to perform the perturbation of priors and likelihood
  v <- c("l_2_p1[1]", "k_2_p1_p2", "k_2_p2_p2", "k_2_p3_p2", "k_2_p4_p2", "phi_2_p2_p1_p2", "phi_2_p2_p2_p2", "phi_2_p2_p3_p2", "phi_2_p2_p4_p2", "alpha_2_p1[1]", "l_2_p1_p2", "l_2_p2_p2", "l_2_p3_p2", "l_2_p4_p2")

  # plot of sensitivity of quantities
  plot_file <- paste0("figures/leaching_priorsense_1_", id_fit, ".pdf")

  priorsense_sequence <-
    priorsense::powerscale_sequence(
      x_stan_fit,
      variable = v,
      moment_match = FALSE
    )

  ggplot2::ggsave(
    plot_file,
    plot =
      priorsense::powerscale_plot_quantities(
        priorsense_sequence,
        variables = v
      ) +
      ggplot2::theme(
        strip.background.x = ggplot2::element_blank(),
        strip.background.y = ggplot2::element_blank()
      ),
    width = width, height = height, dpi = 300,
    device = cairo_pdf,
    limitsize = FALSE
  )


  # get sensivitiy indices
  list(
    powerscale_sensitivity =
      priorsense::powerscale_sensitivity(
        x_stan_fit,
        variable = v,
        lower_alpha = 0.99,
        upper_alpha = 1.01,
        div_measure = "cjs_dist",
        component = c("prior", "likelihood"),
        sensitivity_threshold = 0.05,
        is_method = "psis",
        moment_match = FALSE,
        k_threshold = 0.5,
        resample = FALSE
      ),
    plot_quantities = plot_file
  )

}
