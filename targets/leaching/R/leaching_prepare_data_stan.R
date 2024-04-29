#' Helper function to filter data from Stan and define factors where necessary
#'
#' @param x The main data table
#'
#' @param ...
leaching_prepare_data_stan_1_1 <- function(x, ...) {

  .dots <- list(...)

  res <- x

  # filter
  if(! .dots$has_bengtsson2017) {
    res <-
      res %>%
      dplyr::filter(id_citation != "Bengtsson.2017")
  }

  # define indices
  res <-
    res %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("id_sample_incubation_start", "taxon_rank_value", "speciesxstudies")),
        function(.x) {
          factor(.x, levels = unique(.x))
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
leaching_prepare_data_stan <- function(x, ..., mass_relative_mass_offset) {

  .dots <- list(...)

  res <- leaching_prepare_data_stan_1_1(x = x, ...)

  res <-
    tibble::lst(
      # controls
      has_leaching = .dots$has_leaching,
      has_alpha = .dots$has_alpha,
      prior_only =  .dots$prior_only,
      # data
      N = nrow(res),
      N_has_precision = sum(! is.na(res$mass_relative_mass_precision)),
      m = res$mass_relative_mass,
      t = res$incubation_duration/365,
      precision =
        res %>%
        dplyr::filter(! is.na(mass_relative_mass_precision)) %>%
        dplyr::pull(mass_relative_mass_precision),
      # indices
      index_measurement_has_precision =
        which(! is.na(res$mass_relative_mass_precision)),
      N_samples =
        res %>%
        dplyr::pull(id_sample_incubation_start) %>%
        unique() %>%
        length(),
      index_samples_to_measurements =
        res %>%
        dplyr::pull(id_sample_incubation_start),
      N_species =
        res %>%
        dplyr::pull(taxon_rank_value) %>%
        unique() %>%
        length(),
      index_species_to_samples =
        res %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
        dplyr::pull(taxon_rank_value),
      N_speciesxstudies =
        res %>%
        dplyr::pull(speciesxstudies) %>%
        unique() %>%
        length(),
      index_speciesxstudies_to_samples =
        res %>%
        dplyr::filter(! duplicated(id_sample_incubation_start)) %>%
        dplyr::pull(speciesxstudies),
      # parameters
      k_2_p1_p1 = -2.9,
      k_2_p1_p2 = 0.5,
      phi_2_p2_p1_p1 = 5,
      phi_2_p2_p1_p2 = 0.4,
      k_2_p2_p1 = 0.0,
      k_2_p2_p2 = 0.3,
      phi_2_p2_p2_p1 = 0,
      phi_2_p2_p2_p2 = 0.5,
      k_2_p3_p1 = 0.0,
      k_2_p3_p2 = 0.3,
      phi_2_p2_p3_p1 = 0,
      phi_2_p2_p3_p2 = 0.5,
      phi_2_p2_p4_p1 = 0,
      phi_2_p2_p3_p2 = 0.5,
      k_2_p4_p1 = 0,
      k_2_p4_p2 = 0.3,
      phi_2_p3 = 20,
      phi_2_p1_p1 = 2,
      phi_2_p1_p2 = 2/10,
      phi_2_p1_p3 = 10,
      # leaching
      l_2_p1_p1 = -3.5,
      l_2_p1_p2 = 0.3,
      l_2_p2_p1 = 0.0,
      l_2_p2_p2 = 0.4,
      l_2_p3_p1 = 0.0,
      l_2_p3_p2 = 0.4,
      l_2_p4_p1 = 0,
      l_2_p4_p2 = 0.4,
      # alpha
      alpha_2_p1_p1 = -0.2, # 0.0,
      alpha_2_p1_p2 = 0.3,
      alpha_2_p2_p1 = 0.0,
      alpha_2_p2_p2 = 0.3,
      alpha_2_p3_p1 = 0.0,
      alpha_2_p3_p2 = 0.3,
      alpha_2_p4_p1 = 0,
      alpha_2_p4_p2 = 0.2,
      # group-level standard deviations
      k_2_p1_p2_p1 = 0.4, # 0.3,
      k_2_p2_p2_p1 = 0.4, # 0.3,
      k_2_p3_p2_p1 = 0.4, # 0.3,
      k_2_p4_p2_p1 = 0.4, # 0.3,
      phi_2_p2_p1_p2_p1 = 0.3, # 0.3,
      phi_2_p2_p2_p2_p1 = 0.3, # 0.3,
      phi_2_p2_p3_p2_p1 = 0.3, # 0.3,
      phi_2_p2_p4_p2_p1 = 0.3, # 0.3,
      l_2_p1_p2_p1 = 0.4, # 0.4,
      l_2_p2_p2_p1 = 0.4, # 0.4,
      l_2_p3_p2_p1 = 0.4, # 0.4,
      l_2_p4_p2_p1 = 0.4, # 0.4,
      # alpha_2_p2_p2_p1 = 0.2,
      # alpha_2_p3_p2_p1 = 0.2,
      # constants
      s = mass_relative_mass_offset
    )

}
