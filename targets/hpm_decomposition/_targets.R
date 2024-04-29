# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(magrittr)
library(ggplot2)
library(posterior)

# Set target options:
tar_option_set(
  packages = c("tibble"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# MCMC settings
leaching_mcmc_config <-
  tibble::lst(
    niter = 4000,
    nwarmup = niter/2,
    nchains = 4
  )

# plotting options
ggplot2::theme_set(ggplot2::theme_classic())

#### plan models ####

## real data

# posterior
hpmd_d_models <-
  tibble::tibble(
    data_origin = rlang::syms(c("leaching_stan_draws_4", rep("hpmd_data", 4))),
    has_bengtsson2017 = FALSE,
    has_leaching = TRUE,
    has_alpha = TRUE,
    has_combined_models = c(FALSE, rep(TRUE, 4)),
    has_peat_properties = c(rep(FALSE, 1), rep(TRUE, 4)),
    has_hpm_parameters = c(rep(FALSE, 2), rep(TRUE, 3)),
    has_hpm_leaching = c(rep(FALSE, 3), rep(TRUE, 2)),
    has_cross_validation = c(rep(FALSE, 4), TRUE),
    # has_outlier = TRUE,
    prior_only = c(FALSE)
  ) %>%
  dplyr::mutate(
    id_fit = seq_len(nrow(.)),
    model_name =
      dplyr::case_when(
        id_fit != 5 ~ paste0("3-", id_fit),
        TRUE ~ "3-4"
      )
  )

# priors
hpmd_d_models <-
  dplyr::bind_rows(
    hpmd_d_models,
    hpmd_d_models %>%
      dplyr::slice(2:4) %>%
      dplyr::mutate(
        id_fit = seq_len(nrow(.)) + max(hpmd_d_models$id_fit),
        prior_only = TRUE
      )
  )

# id_fit 4 without outliers
hpmd_d_models <-
  dplyr::bind_rows(
    hpmd_d_models,
    hpmd_d_models %>%
      dplyr::slice(4) %>%
      dplyr::mutate(
        id_fit = seq_len(nrow(.)) + max(hpmd_d_models$id_fit),
        # has_outlier = FALSE,
        model_name = "3-5"
      )
  )

hpmd_has_outlier <- rep(TRUE, nrow(hpmd_d_models))
hpmd_has_outlier[hpmd_d_models$id_fit %in% c(9)] <- FALSE

mass_relative_mass_offset <- 0.0001


# R scripts
lapply(list.files("../leaching/R", full.names = TRUE, recursive = TRUE), source)
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

leaching_d_models <-
  hpmd_d_models |>
  dplyr::mutate(
    model_name =
      model_name |>
      hpmd_model_id_to_name() |>
      factor(levels = paste0(3, "-", 1:5) |> hpmd_model_id_to_name())
  )


# target list
list(
  tar_target(
    leaching_data_from_database_file,
    command = "../leaching/_targets/objects/leaching_data_from_database",
    format = "file"
  ),
  tar_target(
    leaching_data_from_database,
    command = readRDS(leaching_data_from_database_file)
  ),
  tar_target(
    leaching_data_sphagnum_niches_file,
    command = "../leaching/_targets/objects/leaching_data_sphagnum_niches",
    format = "file"
  ),
  tar_target(
    leaching_data_sphagnum_niches,
    command = readRDS(leaching_data_sphagnum_niches_file)
  ),
  tar_target(
    hpmd_data_hpm_microhabitat,
    command =
      hpmd_define_hpm_microhabitats(
        leaching_data_sphagnum_niches = leaching_data_sphagnum_niches
      )
  ),
  tar_target(
    hpmd_data,
    command =
      hpmd_prepare_data(
        x = leaching_data_from_database,
        hpmd_data_hpm_microhabitat = hpmd_data_hpm_microhabitat
      )
  ),
  tar_target(
    leaching_stan_fit_4,
    command = "../leaching/_targets/objects/leaching_stan_fit_4",
    format = "file"
  ),
  tar_target(
    leaching_stan_draws_4,
    command = "../leaching/_targets_rvars/leaching_stan_draws_4.rds",
    format = "file"
  ),
  tar_target(
    hpmd_stan_draws_1,
    command =
      hpmd_get_yhat_fit_1(
        x = readRDS_rvars(leaching_stan_draws_4),
        hpmd_data_hpm_microhabitat = hpmd_data_hpm_microhabitat,
        id_fit = hpmd_d_models$id_fit[[1]]
      ),
    format = "file"
  ),
  tar_target(
    hpmd_stan_rmse_1,
    command = {
      cur_d_models <- hpmd_d_models[hpmd_d_models$id_fit == 1, ]
      res <- paste0("_targets_rvars/hpmd_stan_rmse_train_", 1, ".rds")
      hpmd_get_rmse(
        x_stan_draws = readRDS_rvars(hpmd_stan_draws_1),
        has_bengtsson2017 = cur_d_models$has_bengtsson2017[[1]],
        has_leaching = cur_d_models$has_leaching[[1]],
        has_alpha = cur_d_models$has_alpha[[1]],
        has_hpm_parameters = cur_d_models$has_hpm_parameters[[1]],
        has_hpm_leaching = cur_d_models$has_hpm_leaching[[1]],
        has_cross_validation = cur_d_models$has_cross_validation[[1]],
        mass_relative_mass_offset = mass_relative_mass_offset
      ) %>%
        saveRDS_rvars(file = res)
      res
    },
    format = "file"
  ),
  tar_target(
    hpmd_mm39_file,
    command = {
      "stan/mm39.stan"
    },
    format = "file"
  ),
  tar_target(
    hpmd_mm39,
    command = {
      hpmd_mm39_file
      rstan::stan_model("stan/mm39.stan", model_name = "mm39")
    }
  ),
  tarchetypes::tar_map(
    values = list(cur_id_fit = hpmd_d_models$id_fit[-c(1)], cur_data = hpmd_d_models$data_origin[-c(1)]),
    names = dplyr::all_of("cur_id_fit"),
    tar_target(
      hpmd_stan_cur_d_models,
      command = hpmd_d_models[hpmd_d_models$id_fit == cur_id_fit, ]
    ),
    tar_target(
      hpmd_data_stan,
      command = {
        hpmd_prepare_data_stan(
          x = eval(cur_data),
          has_bengtsson2017 = hpmd_stan_cur_d_models$has_bengtsson2017[[1]],
          has_leaching = hpmd_stan_cur_d_models$has_leaching[[1]],
          has_alpha = hpmd_stan_cur_d_models$has_alpha[[1]],
          prior_only =  hpmd_stan_cur_d_models$prior_only[[1]],
          has_hpm_parameters = hpmd_stan_cur_d_models$has_hpm_parameters[[1]],
          has_hpm_leaching = hpmd_stan_cur_d_models$has_hpm_leaching[[1]],
          has_cross_validation = hpmd_stan_cur_d_models$has_cross_validation[[1]],
          has_outlier = hpmd_has_outlier[[cur_id_fit]],
          mass_relative_mass_offset = mass_relative_mass_offset,
          hpmd_stan_draws_1 = hpmd_stan_draws_1
        )
      }
    ),
    tar_target(
      hpmd_stan_initial_values,
      command = {
        purrr::map(seq_len(leaching_mcmc_config$nchains), function(i) {
          res <-
            list(
              hpm_k_2_p1 = runif(1, 0.8, 1.5),
              phi_2_p2_p1_p2 = runif(1, 0, 4),
              k_2_p1_p2 = runif(1, 0, 2),
              k_2_p1 = runif(1, 0, 2),
              phi_2_p1 = runif(1, 0, 1)
            )
          if(hpmd_stan_cur_d_models$has_hpm_parameters[[1]]) {
            res$m68_p2 = array(runif(1, 0, 1.5), dim = 1L)
          }
          if(hpmd_stan_cur_d_models$has_hpm_leaching[[1]]) {
            res$hpm_l_2_p4 = array(runif(1, 30, 50), dim = 1L)
          }
          res
        })
      }
    ),
    tar_target(
      hpmd_stan_fit,
      command = {
        rstan::sampling(
          hpmd_mm39,
          data = tidybayes::compose_data(hpmd_data_stan),
          iter = leaching_mcmc_config$niter,
          warmup = leaching_mcmc_config$nwarmup,
          chains = leaching_mcmc_config$nchains,
          cores = leaching_mcmc_config$nchains,
          algorithm = "NUTS",
          control = list(max_treedepth = 14, adapt_delta = 0.99),
          init = hpmd_stan_initial_values
        )
      }
    ),
    tar_target(
      hpmd_stan_fit_cv,
      command = {
        if(hpmd_stan_cur_d_models$has_cross_validation[[1]]) {
          purrr::map(seq_along(hpmd_data_stan), function(i) {
            print(i)
            rstan::sampling(
              hpmd_mm39,
              data = tidybayes::compose_data(hpmd_data_stan[[i]]),
              iter = leaching_mcmc_config$niter,
              warmup = leaching_mcmc_config$nwarmup,
              chains = leaching_mcmc_config$nchains,
              cores = leaching_mcmc_config$nchains,
              algorithm = "NUTS",
              control = list(max_treedepth = 14, adapt_delta = 0.99),
              init = hpmd_stan_initial_values
            )
          })
        } else {
          NULL
        }
      }
    ),
    tar_target(
      hpmd_stan_draws,
      command = {
        res <- paste0("_targets_rvars/hpmd_stan_draws_", cur_id_fit, ".rds")
        hpmd_stan_extract_draws(
          x_stan_model =
            if(hpmd_stan_cur_d_models$has_cross_validation[[1]]) {
              hpmd_stan_fit_cv
            } else {
              hpmd_stan_fit
            },
          x = eval(cur_data),
          x_stan_data = hpmd_data_stan,
          data_template =
            eval(cur_data),#  %>%
          #dplyr::select(id_sample, id_sample_incubation_start, taxon_rank_value, speciesxstudies),
          id_fit = cur_id_fit,
          has_bengtsson2017 = hpmd_stan_cur_d_models$has_bengtsson2017[[1]],
          has_leaching = hpmd_stan_cur_d_models$has_leaching[[1]],
          has_alpha = hpmd_stan_cur_d_models$has_alpha[[1]],
          has_hpm_parameters = hpmd_stan_cur_d_models$has_hpm_parameters[[1]],
          has_hpm_leaching = hpmd_stan_cur_d_models$has_hpm_leaching[[1]],
          has_cross_validation = hpmd_stan_cur_d_models$has_cross_validation[[1]],
          has_outlier = hpmd_has_outlier[[cur_id_fit]],
          mass_relative_mass_offset = mass_relative_mass_offset,
          hpmd_stan_draws_1 = hpmd_stan_draws_1
        ) %>%
          saveRDS_rvars(file = res)
        res
      },
      format = "file"
    ),
    tar_target(
      hpmd_stan_rmse,
      command = {
        res <- paste0("_targets_rvars/hpmd_stan_rmse_train_", cur_id_fit, ".rds")
        hpmd_get_rmse(
          x_stan_draws = readRDS_rvars(hpmd_stan_draws),
          has_bengtsson2017 = hpmd_stan_cur_d_models$has_bengtsson2017[[1]],
          has_leaching = hpmd_stan_cur_d_models$has_leaching[[1]],
          has_alpha = hpmd_stan_cur_d_models$has_alpha[[1]],
          has_hpm_parameters = hpmd_stan_cur_d_models$has_hpm_parameters[[1]],
          has_hpm_leaching = hpmd_stan_cur_d_models$has_hpm_leaching[[1]],
          has_cross_validation = hpmd_stan_cur_d_models$has_cross_validation[[1]],
          mass_relative_mass_offset = mass_relative_mass_offset
        ) %>%
          saveRDS_rvars(file = res)
        res
      },
      format = "file"
    ),
    tar_target(
      hpmd_plot_4,
      command =
        if(cur_id_fit %in% c(3:4, 9)) {
          hpmd_get_plot_4(
            x_stan_draws = readRDS_rvars(hpmd_stan_draws)
          )
        } else {
          NA_character_
        }
    )
  ),
  tar_target(
    hpmd_stan_checks,
    command =
      hpmd_stan_get_all_checks(
        prior_models =
          list(
            hpmd_stan_draws_6,
            hpmd_stan_draws_7,
            hpmd_stan_draws_8
          ),
        posterior_models =
          list(
            hpmd_stan_draws_1,
            hpmd_stan_draws_2,
            hpmd_stan_draws_3,
            hpmd_stan_draws_4,
            hpmd_stan_draws_9
          ),
        x_stan_models =
          list(
            readRDS(leaching_stan_fit_4),
            hpmd_stan_fit_2,
            hpmd_stan_fit_3,
            hpmd_stan_fit_4,
            hpmd_stan_fit_9
          )
      )
  ),
  tar_target(
    hpmd_plot_1,
    command =
      hpmd_get_plot_1(
        x_stan_draws =
          list(
            hpmd_stan_draws_1,
            hpmd_stan_draws_2,
            hpmd_stan_draws_3,
            hpmd_stan_draws_4
          )
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_2,
    command =
      hpmd_get_plot_2(
        x_stan_draws =
          list(
            hpmd_stan_draws_1,
            hpmd_stan_draws_2,
            hpmd_stan_draws_3,
            hpmd_stan_draws_4
          )
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_3,
    command =
      hpmd_get_plot_3(
        x_stan_draws =
          list(
            hpmd_stan_draws_1,
            hpmd_stan_draws_2,
            hpmd_stan_draws_3,
            hpmd_stan_draws_4
          )
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_5,
    command =
      hpmd_get_plot_5(
        x_stan_draws = readRDS_rvars(hpmd_stan_draws_5)
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_6,
    command =
      hpmd_get_plot_6(
        x_stan_draws =
          list(
            readRDS_rvars(hpmd_stan_draws_3),
            readRDS_rvars(hpmd_stan_draws_4)
          ),
        hpmd_d_models = hpmd_d_models
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_7,
    command =
      hpmd_get_plot_7(
        x_stan_draws =
          list(
            hpmd_stan_draws_1,
            hpmd_stan_draws_2,
            hpmd_stan_draws_3,
            hpmd_stan_draws_4
          )
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_8_1,
    command =
      hpmd_get_plot_8(
        x_stan_draws = readRDS_rvars(hpmd_stan_draws_1),
        file_plot = "figures/hpmd_plot_8_1.pdf"
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_8_4,
    command =
      hpmd_get_plot_8(
        x_stan_draws = readRDS_rvars(hpmd_stan_draws_4),
        file_plot = "figures/hpmd_plot_8_4.pdf"
      ),
    format = "file"
  ),
  tar_target(
    hpmd_plot_9,
    command =
      hpmd_get_plot_9(
        x_stan_draws =
          list(
            hpmd_stan_draws_2,
            hpmd_stan_draws_3,
            hpmd_stan_draws_4
          ),
        hpmd_stan_draws_1 = readRDS_rvars(hpmd_stan_draws_1)
      ),
    format = "file"
  ),
  tar_map(
    values = list(cur_id_simulation = 1:6),
    names = dplyr::all_of("cur_id_simulation"),
    tar_target(
      hpmd_simulation_1_newdata,
      command =
        hpmd_get_newdata_for_simulations(
          id_simulation = cur_id_simulation,
          x_stan_fit_4 = hpmd_stan_fit_4,
          x_stan_data_4 = hpmd_data_stan_4,
          hpmd_data = hpmd_data
        ),
      format = "file"
    ),
    tar_target(
      hpmd_simulation_1_plot,
      command =
        hpmd_get_plots_for_simulations(
          newdata = readRDS_rvars(hpmd_simulation_1_newdata),
          id_simulation = cur_id_simulation,
          x_stan_fit_4 = hpmd_stan_fit_4,
          x_stan_data_4 = hpmd_data_stan_4,
          hpmd_data = hpmd_data
        ),
      format = "file"
    )
  ),
  tar_target(
    hpmd_tables,
    command =
      list(
        prior_justification =
          hpmd_get_table_prior_justification(
            x_stan_data = hpmd_data_stan_4,
            nsim = with(leaching_mcmc_config, (niter - nwarmup) * nchains)
          )
      )
  ),
  tar_target(
    hpmd_priorsense_fit_4,
    command = {
      v <- c("l_2_p1", "k_2_p1", "alpha_2_p1", "phi_2_p1", "layer_total_porosity_1", "layer_minimum_degree_of_saturation_at_surface_1", "layer_water_table_depth_to_surface_1", "hpm_k_2_p1", "m69_p1", "m69_p2", "m68_p1", "m68_p2", "m68_p3_2_p1", "hpm_l_2_p4", "hpm_l_2_p1", "hpm_l_2_p3")
      priorsense::powerscale_sensitivity(hpmd_stan_fit_4, variable = v)
    }
  ),
  tar_render(
    hpmd_supporting_info,
    path = "hpmd-supporting-info.Rmd"
  ),
  tar_render(
    hpmd_paper,
    path = "hpmd-paper.Rmd"
  )
)




