# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(magrittr)
library(ggplot2)
library(priorsense)

# Set target options:
tar_option_set(
  packages = c("RMariaDB", "dpeatdecomposition", "tibble", "dplyr", "magrittr", "rstan", "ggplot2", "ggdist", "posterior", "dbplyr", "purrr", "dm", "ggplot2"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

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
leaching_d_models <-
  tidyr::expand_grid(
    data_origin = rlang::syms("leaching_data_from_database"),
    has_leaching = c(TRUE, FALSE),
    has_bengtsson2017 = c(TRUE, FALSE),
    has_alpha = c(TRUE, FALSE),
    prior_only = c(FALSE)
  ) %>%
  dplyr::arrange(dplyr::desc(has_leaching), has_alpha, dplyr::desc(has_bengtsson2017)) %>%
  dplyr::mutate(
    id_fit = seq_len(nrow(.))
  ) %>%
  dplyr::group_by(has_leaching) %>%
  dplyr::mutate(
    model_name =
      paste0(
        as.integer(1 - has_leaching + 1L), "-",
        # as.integer(as.factor(has_bengtsson2017)), "-",
        # as.integer(as.factor(has_alpha))
        dplyr::consecutive_id(seq_len(dplyr::n()))
      )
  ) %>%
  dplyr::ungroup()

# prior
leaching_d_models <-
  dplyr::bind_rows(
    leaching_d_models,
    leaching_d_models %>%
      dplyr::filter(! has_bengtsson2017) %>%
      dplyr::filter(! duplicated(paste0(has_leaching, "_", has_alpha))) %>%
      dplyr::mutate(
        prior_only = TRUE,
        id_fit = seq_len(nrow(.)) + max(leaching_d_models$id_fit),
        model_name =
          paste0(
            as.integer(1 - has_leaching + 1L), "-",
            # as.integer(as.factor(has_bengtsson2017)), "-",
            # as.integer(as.factor(has_alpha))
            rep(c(2, 4), 2)
          )
      )
  )

## simulated data

# leaching_data_sim_1
leaching_d_models <-
  dplyr::bind_rows(
    leaching_d_models,
    tibble::tibble(
      data_origin = rlang::syms("leaching_data_sim_1"),
      has_leaching = c(TRUE),
      has_bengtsson2017 = c(FALSE),
      has_alpha = c(TRUE),
      prior_only = c(FALSE),
      id_fit = seq_along(prior_only) + max(leaching_d_models$id_fit),
      model_name =
        paste0(
          as.integer(1 - has_leaching + 1L), "-",
          # as.integer(as.factor(has_bengtsson2017)), "-",
          # as.integer(as.factor(has_alpha))
          5
        )
    )
  )


# leaching_data_sim_2
leaching_d_models <-
  dplyr::bind_rows(
    leaching_d_models,
    tibble::tibble(
      data_origin = rlang::syms("leaching_data_sim_2"),
      has_leaching = c(TRUE),
      has_bengtsson2017 = c(FALSE),
      has_alpha = c(TRUE),
      prior_only = c(FALSE),
      id_fit = seq_along(prior_only) + max(leaching_d_models$id_fit),
      model_name =
        paste0(
          as.integer(1 - has_leaching + 1L), "-",
          # as.integer(as.factor(has_bengtsson2017)), "-",
          # as.integer(as.factor(has_alpha))
          6
        )
    )
  )

mass_relative_mass_offset <- 0.0001

# leaching_test_1 <- c(1, 2)

# Load the R scripts with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)


# Target list
list(
  tar_target(
    leaching_data_sphagnum_niches,
    command =
      leaching_data_get_niche_wtd(file = "derived_data/evo12547-sup-0001-suppmat(1)-Tab-S2.csv")
  ),
  tar_target(
    leaching_dpeatdecomposition_snapshot,
    command = "derived_data/dpeatdecomposition_snapshot.rds",
    format = "file"
  ),
  tar_target(
    leaching_data_from_database,
    command =
      leaching_prepare_data(
        file = leaching_dpeatdecomposition_snapshot,
        mass_relative_mass_offset = mass_relative_mass_offset,
        leaching_data_sphagnum_niches = leaching_data_sphagnum_niches
      )
  ),
  tar_target(
    leaching_data_sim_1,
    command =
      leaching_get_simulated_data(
        id_dataset = "leaching_data_sim_1",
        mass_relative_mass_offset = mass_relative_mass_offset
      )
  ),
  tar_target(
    leaching_data_sim_2,
    command =
      leaching_get_simulated_data(
        id_dataset = "leaching_data_sim_2",
        mass_relative_mass_offset = mass_relative_mass_offset,
        x = leaching_data_from_database,
        x_model_draws = readRDS_rvars(leaching_stan_draws_4)
      )
  ),
  tar_target(
    leaching_mm38_file,
    command = {
      "stan/mm38.stan"
    },
    format = "file"
  ),
  tar_target(
    leaching_mm38,
    command = {
      leaching_mm38_file
      rstan::stan_model("stan/mm38.stan", model_name = "mm38")
    }
  ),
  tarchetypes::tar_map(
     values = list(cur_id_fit = leaching_d_models$id_fit, cur_data = leaching_d_models$data_origin),
     names = dplyr::all_of("cur_id_fit"),
     tar_target(
       leaching_data_stan,
       command = {
           cur_d_models <- leaching_d_models[leaching_d_models$id_fit == cur_id_fit, ]
           leaching_prepare_data_stan(
             x = eval(cur_data),
             has_bengtsson2017 = cur_d_models$has_bengtsson2017[[1]],
             has_leaching = cur_d_models$has_leaching[[1]],
             has_alpha = cur_d_models$has_alpha[[1]],
             prior_only =  cur_d_models$prior_only[[1]],
             mass_relative_mass_offset = mass_relative_mass_offset
           )
         }
     ),
     tar_target(
       leaching_stan_fit,
       command = {
         rstan::sampling(
           leaching_mm38,
           data = tidybayes::compose_data(leaching_data_stan),
           iter = leaching_mcmc_config$niter,
           warmup = leaching_mcmc_config$nwarmup,
           chains = leaching_mcmc_config$nchains,
           cores = leaching_mcmc_config$nchains,
           algorithm = "NUTS",
           control = list(max_treedepth = 12, adapt_delta = 0.99)
         )
       }
     ),
     tar_target(
       leaching_stan_draws,
       command = {
         cur_d_models <- leaching_d_models[leaching_d_models$id_fit == cur_id_fit, ]
         res <- paste0("_targets_rvars/leaching_stan_draws_", cur_id_fit, ".rds")
         leaching_stan_extract_draws(
           x_stan_model = leaching_stan_fit,
           x = eval(cur_data),
           x_stan_data = leaching_data_stan,
           data_template =
             eval(cur_data),#  %>%
             #dplyr::select(id_sample, id_sample_incubation_start, taxon_rank_value, speciesxstudies),
           id_fit = cur_id_fit,
           has_bengtsson2017 = cur_d_models$has_bengtsson2017[[1]],
           has_leaching = cur_d_models$has_leaching[[1]],
           has_alpha = cur_d_models$has_alpha[[1]]
         ) %>%
           saveRDS_rvars(file = res)
         res
       },
       format = "file"
     ),
     tar_target(
       leaching_plot_1,
       command = {
         cur_d_models <- leaching_d_models[leaching_d_models$id_fit == cur_id_fit, ]
         leaching_plot_estimates_for_species_and_studies(
           x_stan_draws = readRDS_rvars(leaching_stan_draws),
           id_fit = cur_id_fit,
           has_leaching = cur_d_models$has_leaching[[1]],
           has_alpha = cur_d_models$has_alpha[[1]]
         )
       }
     ),
     tar_target(
       leaching_uncertainty_analysis_1,
       command = {
         cur_d_models <- leaching_d_models[leaching_d_models$id_fit == cur_id_fit, ]
         leaching_make_uncertainty_analysis_1(
           x_stan_draws = readRDS_rvars(leaching_stan_draws),
           nbin = 10,
           data_template = eval(cur_data),
           has_leaching = cur_d_models$has_leaching[[1]],
           has_alpha = cur_d_models$has_alpha[[1]]
         )
       }
     )
  ),
  tar_target(
    leaching_stan_checks,
    command =
      leaching_stan_get_all_checks(
        prior_models =
          list(
            leaching_stan_draws_9,
            leaching_stan_draws_10,
            leaching_stan_draws_11,
            leaching_stan_draws_12
          ),
        posterior_models =
          list(
            leaching_stan_draws_1,
            leaching_stan_draws_2,
            leaching_stan_draws_3,
            leaching_stan_draws_4,
            leaching_stan_draws_5,
            leaching_stan_draws_6,
            leaching_stan_draws_7,
            leaching_stan_draws_8,
            leaching_stan_draws_13,
            leaching_stan_draws_14
          ),
        x_stan_models =
          list(
            leaching_stan_fit_1,
            leaching_stan_fit_2,
            leaching_stan_fit_3,
            leaching_stan_fit_4,
            leaching_stan_fit_5,
            leaching_stan_fit_6,
            leaching_stan_fit_7,
            leaching_stan_fit_8,
            leaching_stan_fit_13,
            leaching_stan_fit_14
          )
      )
  ),
  tar_target(
    leaching_plot_2_fit_5_and_fit_1,
    command =
      leaching_plot_k0_l0(
        x_stan_draws_no_leaching = readRDS_rvars(leaching_stan_draws_5),
        x_stan_draws_with_leaching = readRDS_rvars(leaching_stan_draws_1),
        n_experiments_per_species = 5,
        facet_nrow = 3,
        width = 8,
        height = 9
      ),
    format = "file"
  ),
  tar_target(
    leaching_plot_2_fit_6_and_fit_2,
    command =
      leaching_plot_k0_l0(
        x_stan_draws_no_leaching = readRDS_rvars(leaching_stan_draws_6),
        x_stan_draws_with_leaching = readRDS_rvars(leaching_stan_draws_2),
        n_experiments_per_species = 5,
        facet_nrow = 3,
        width = 8,
        height = 9
      ),
    format = "file"
  ),
  tar_target(
    leaching_plot_2_fit_7_and_fit_3,
    command =
      leaching_plot_k0_l0(
        x_stan_draws_no_leaching = readRDS_rvars(leaching_stan_draws_7),
        x_stan_draws_with_leaching = readRDS_rvars(leaching_stan_draws_3),
        n_experiments_per_species = 5,
        facet_nrow = 3,
        width = 8,
        height = 9
      ),
    format = "file"
  ),
  tar_target(
    leaching_plot_3,
    command =
      leaching_get_plot_3(
        x_ua_1 =
          leaching_uncertainty_analysis_1_4 |>
          dplyr::filter(id_citation != "Bengtsson.2017"),
        x_stan_draws =
          readRDS_rvars(leaching_stan_draws_4) |>
          dplyr::filter(id_citation != "Bengtsson.2017")
      ),
    format = "file"
  ),
  tar_target(
    leaching_plot_4,
    command =
      list(
        fit_1_vs_fit_2 =
          leaching_get_plot_4_individual(
            x_stan_draws = readRDS_rvars(leaching_stan_draws_1),
            y_stan_draws = readRDS_rvars(leaching_stan_draws_2)
          ),
        fit_3_vs_fit_4 =
          leaching_get_plot_4_individual(
            x_stan_draws = readRDS_rvars(leaching_stan_draws_3),
            y_stan_draws = readRDS_rvars(leaching_stan_draws_4)
          )
      )
  ),
  tar_target(
    leaching_plot_5,
    command = leaching_get_plot_5(readRDS_rvars(leaching_stan_draws_13))
  ),
  tar_target(
    leaching_plot_outlier_1_3,
    command =
      leaching_get_plot_outlier_trajectories(
        x_stan_draws = readRDS_rvars(leaching_stan_draws_3),
        alpha = 0.99,
        width = 10,
        height = 8
        )
  ),
  tar_target(
    leaching_simulation_interpretation_alpha_1,
    command =
      leaching_get_simulation_interpretation_alpha_1()
  ),
  tar_target(
    leaching_tables,
    command =
      list(
        prior_justification =
          leaching_get_table_prior_justification(
            x_stan_data = leaching_data_stan_1,
            nsim = with(leaching_mcmc_config, (niter - nwarmup) * nchains)
          ),
        parameter_estimates =
          leaching_get_table_estimates(
            x_stan_draws =
              list(
                leaching_stan_draws_1,
                leaching_stan_draws_2,
                leaching_stan_draws_3,
                leaching_stan_draws_4,
                leaching_stan_draws_5,
                leaching_stan_draws_6,
                leaching_stan_draws_7,
                leaching_stan_draws_8
              )
          ),
        parameter_estimates_species =
          leaching_get_table_estimates_species(
            x_stan_draws =
              list(
                leaching_stan_draws_3,
                leaching_stan_draws_4
              )
          )
      )
  ),
  tar_target(
    leaching_priorsense_1_4,
    command =
      leaching_do_priorsense_analysis(
        x_stan_fit = leaching_stan_fit_4,
        id_fit = 4,
        width = 30,
        height = 50
      )
  ),
  tar_render(
    leaching_paper,
    path = "leaching-paper.Rmd"
  ),
  tar_render(
    leaching_supporting_info,
    path = "leaching-supporting-info.Rmd"
  )
)
