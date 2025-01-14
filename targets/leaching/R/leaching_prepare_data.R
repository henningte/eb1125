#### Main table ####

#' Extracts data from the Sphagnum Decomposition Database
#'
#' @param ... Arguments passed to `RMariaDB::dbConnect()`.
#'
#' @export
leaching_get_data_database <- function(..., file) {


  # connect to database
  con <-
    RMariaDB::dbConnect(
      drv = RMariaDB::MariaDB(),
      dbname = "dpeatdecomposition",
      ...
    )

  # get database as dm object
  dm_dpeatdecomposition <-
    dpeatdecomposition::dp_get_dm(con, learn_keys = TRUE)

  # get information on samples and measurements
  res <-
    dm_dpeatdecomposition %>%
    dm::dm_zoom_to(samples) %>%
    dm::left_join(data, by = "id_sample") %>%
    dm::pull_tbl() %>%
    tibble::as_tibble()

  citations_to_datasets <-
    dm_dpeatdecomposition %>%
    dm::dm_zoom_to(citations_to_datasets) %>%
    dm::pull_tbl() %>%
    tibble::as_tibble()

  # pmird_units <-
  #   dm_dpeatdecomposition %>%
  #   dm::dm_zoom_to(attributes) %>%
  #   dm::left_join(measurement_scales, by = "id_measurement_scale") %>%
  #   dm::left_join(measurement_scales_ratio, by = "id_measurement_scale") %>%
  #   dm::left_join(units, by = "id_unit") %>%
  #   dm::pull_tbl() %>%
  #   tibble::as_tibble()

  # disconnect
  on.exit(RMariaDB::dbDisconnect(con), add = TRUE)

  # add bibkeys to identify datasets
  res <-
    dplyr::left_join(
      res,
      citations_to_datasets %>%
        dplyr::filter(! duplicated(id_dataset)),
      by = "id_dataset"
    )

  saveRDS(res, file)

  file

}

#' Loads the data extracted with `leaching_get_data_database()`
#'
#' @param file Character value. Path to the file extracted with `leaching_get_data_database()`.
#'
#' @export
leaching_load_data_database <- function(file) {
  readRDS(file)
}


#' Filters the data
#'
#' @param x Output from `leaching_load_data_database()`.
#'
#' @export
leaching_filter_data_1 <- function(x) {

  # which of the datasets used Sphagnum litter?
  x_has_sphagnum <-
    x %>%
    dplyr::group_by(id_dataset) %>%
    dplyr::summarise(
      dataset_has_sphagnum =
        any(stringr::str_detect(taxon_rank_value, "Sphagnum")),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dataset_has_sphagnum = ifelse(is.na(dataset_has_sphagnum), FALSE, dataset_has_sphagnum)
    )

  # which of the dataset has information on `"water_table_depth"`?
  x_has_water <-
    x %>%
    dplyr::group_by(id_dataset) %>%
    dplyr::summarise(
      dataset_has_water =
        any(attribute_name %in% c("water_table_depth")),
      .groups = "drop"
    )

  # identify datasets
  res <-
    purrr::reduce(
      list(x, x_has_sphagnum, x_has_water),
      dplyr::left_join,
      by = "id_dataset"
    )

  res

}


#' Adds metadata
#'
#' @param x Output from `leaching_filter_data_1()`.
#'
#' @export
leaching_add_metadata_1 <- function(x) {

  # identify origin sample_type and origin sampling_month
  res <-
    x %>%
    dplyr::mutate(
      origin_sample_type =
        purrr::map_chr(seq_len(nrow(.)), function(i) {
          if(! is.na(id_sample_origin[[i]])) {
            sample_type[id_sample_origin[[i]]]
          } else {
            NA_character_
          }
        }),
      origin_sampling_month =
        purrr::map_chr(seq_len(nrow(.)), function(i) {
          if(! is.na(id_sample_origin[[i]])) {
            sampling_month[id_sample_origin[[i]]]
          } else {
            NA_real_
          }
        }),
      origin_sample_depth_upper =
        purrr::map_chr(seq_len(nrow(.)), function(i) {
          if(! is.na(id_sample_origin[[i]])) {
            sample_depth_upper[id_sample_origin[[i]]]
          } else {
            NA_real_
          }
        }),
      origin_sample_depth_lower =
        purrr::map_chr(seq_len(nrow(.)), function(i) {
          if(! is.na(id_sample_origin[[i]])) {
            sample_depth_lower[id_sample_origin[[i]]]
          } else {
            NA_real_
          }
        })
    )

  res

}


#' Converts to wide format (one row per sample)
#'
#' @param
#'
#' @export
leaching_data_make_wide <- function(x) {

  # variables for which we need values
  target_attribute_names <-
    c("mass_relative_mass", "mesh_size_absolute", "water_table_depth")

  # column types for each variable
  target_column_types <- c("value", "error", "value_type", "error_type", "sample_size")

  # filter for variables
  res <-
    x %>%
    dplyr::filter(attribute_name %in% target_attribute_names) %>%
    dplyr::filter(! duplicated(paste0(id_sample, "_", attribute_name)))

  # widen column types
  res <-
    purrr::map(target_column_types, function(.x) {

      # filter
      res <-
        if(.x == "value") {
          res %>%
            dplyr::select(! dplyr::all_of(c("error", "error_type", "value_type", "sample_size", "comments_measurement", "id_measurement", "id_measurement_numerator", "id_measurement_denominator")))
        } else {
          res %>%
            dplyr::select(dplyr::all_of(c("id_sample", .x, "attribute_name")))
        }

      # widen
      res <-
        res %>%
        tidyr::pivot_wider(
          names_from = attribute_name,
          values_from = .x
        )

      # new column names
      if(.x != "value") {
        res <-
          res %>%
          setNames(nm = c("id_sample", paste0(colnames(.)[-1], "_", .x)))
      }

      res

    }) %>%
    purrr::reduce(dplyr::left_join, by = "id_sample")

  res

}


#' Extracts water table depth data and assigns it to litterbag samples
#'
#' @param
#'
#' @export
leaching_get_average_wtd <- function(x) {

  # target variables
  target_variable_names <- c("sampling_year", "sampling_month", "sampling_day", "sample_depth_upper", "sample_depth_lower", "water_table_depth", "water_table_depth_error", "water_table_depth_value_type", "water_table_depth_error_type", "water_table_depth_sample_size")

  # index to filter rows (either there are water table depth data for each litterbag available or the study reports water table depth separately)
  index <- '! (! is.na(water_table_depth) & ! is.na(mass_relative_mass)) & sample_type == "peat"'

  # extract water table depth and collect them in a list column for each id_dataset x experimental_design
  res <-
    x %>%
    dplyr::filter(eval(rlang::parse_expr(index))) %>%
    dplyr::select(dplyr::all_of(c("id_dataset", "experimental_design", target_variable_names))) %>%
    tidyr::nest(
      data_water_table_depth = dplyr::all_of(target_variable_names)
    )

  # convert errors to sd where possible, otherwise set errors to NA
  res <-
    res %>%
    dplyr::mutate(
      data_water_table_depth =
        purrr::map(data_water_table_depth, function(.x) {
          .x %>%
            dplyr::mutate(
              water_table_depth_error =
                dplyr::case_when(
                  water_table_depth_error_type == "se" ~ water_table_depth_error * sqrt(water_table_depth_sample_size),
                  water_table_depth_error_type == "sd" ~ water_table_depth_error,
                  TRUE ~ NA_real_
                )
            )
        })
    )

  # compute average water table depth
  res <-
    res %>%
    dplyr::mutate(
      water_table_depth_average =
        purrr::map_dbl(data_water_table_depth, function(.x) {
          .x %>%
            dplyr::mutate(
              water_table_depth = errors::set_errors(water_table_depth, water_table_depth_error)
            ) %>%
            dplyr::filter(sample_depth_upper == 0 & sample_depth_lower == 0) %>% #---note: only water table depth relative to the surface
            dplyr::pull(water_table_depth) %>%
            mean(na.rm = TRUE)
        }),
      water_table_depth_average_error = errors::errors(water_table_depth_average),
      water_table_depth_average = as.numeric(water_table_depth_average)
    )

  # assign to x
  res <-
    dplyr::left_join(
      x,
      res %>%
        dplyr::select(id_dataset, experimental_design, water_table_depth_average, water_table_depth_average_error),
      by = c("id_dataset", "experimental_design")
    )

  # when there is already a water table depth for each litterbag, use this value
  res <-
    res %>%
    dplyr::mutate(
      water_table_depth_error =
        dplyr::case_when(
          water_table_depth_error_type == "se" ~ water_table_depth_error * sqrt(water_table_depth_sample_size),
          water_table_depth_error_type == "sd" ~ water_table_depth_error,
          TRUE ~ NA_real_
        ),
      water_table_depth_average =
        dplyr::case_when(
          is.na(water_table_depth) ~ water_table_depth_average,
          TRUE ~ water_table_depth + sample_depth_upper #---note: to compute the water table depth relative to the surface
        ),
      water_table_depth_average_error =
        dplyr::case_when(
          is.na(water_table_depth) ~ water_table_depth_average_error,
          TRUE ~ water_table_depth_error
        )
    ) %>%
    dplyr::rename(
      water_table_depth_to_surface_average = "water_table_depth_average",
      water_table_depth_to_surface_average_error = "water_table_depth_average_error"
    )

  res

}

#' Final filter step: Keep only values, where remaining masses are available and for Sphagnum
#'
#' @param
#'
#' @export
leaching_filter_data_2 <- function(x) {

  # filter index
  x_has_more_than_one_measurement <-
    x %>%
    dplyr::filter(incubation_duration > 0.0 & ! is.na(incubation_duration) & ! is.na(id_sample_incubation_start)) %>%
    dplyr::group_by(id_sample_incubation_start) %>%
    dplyr::summarise(
      has_more_than_one_measurement = length(incubation_duration) > 1
    )

  # filter
  res <-
    x %>%
    dplyr::left_join(x_has_more_than_one_measurement, by = "id_sample_incubation_start") %>%
    dplyr::filter(! is.na(mass_relative_mass) & (stringr::str_detect(taxon_rank_value, "Sphagnum") | sample_type == "peat") & has_more_than_one_measurement & incubation_duration > 0.0)

  # sort
  res <-
    res %>%
    dplyr::arrange(id_sample_incubation_start, incubation_duration)

  res

}


#' Adjusts taxon_rank_value
#'
#' @param x
#'
#' @export
leaching_data_adjust_taxon_rank_value <- function(x) {

  x %>%
    dplyr::mutate(
      # make taxon_rank_value more specific where possible
      taxon_rank_value =
        dplyr::case_when(
          id_citation == "Hagemann.2015" & taxon_rank_value == "Sphagnum" ~ "Sphagnum russowii or capillifolium",
          TRUE ~ taxon_rank_value
        ),
      speciesxstudies = paste0(taxon_rank_value, "_", id_dataset)
    )

}





#' Prepares remaining mass values for modeling
#'
#' @param x
#'
#' @export
leaching_prepare_remaining_masses <- function(x, mass_relative_mass_offset) {

  x %>%
    dplyr::mutate(
      # convert errors to sd
      mass_relative_mass_error =
        dplyr::case_when(
          mass_relative_mass_error_type == "se" ~ mass_relative_mass_error * sqrt(mass_relative_mass_sample_size),
          mass_relative_mass_error_type == "sd" ~ mass_relative_mass_error,
          TRUE ~ NA_real_
        ),
      # correct values > 1
      mass_relative_mass =
        ifelse(
          mass_relative_mass >= 1.0,
          1.0 - mass_relative_mass_offset,
          mass_relative_mass
        ),
      # compute precision for beta distribution from sd where possible
      mass_relative_mass_precision =
        purrr::map_dbl(seq_along(mass_relative_mass), function(i) {
          mu <- mass_relative_mass[[i]]
          sigma <- mass_relative_mass_error[[i]]
          if(is.na(mu) || is.na(sigma)) {
            res <- NA_real_
          } else {
            if(sigma <= 0.0) {
              sigma <- 0.001 #---note: A smaller value causes way too large precision values which cause divergent transitions
            }
            if(sigma^2 >= (mu * (1 - mu))) {
              sigma <- sqrt(mu * (1 - mu)) * (1 - 0.001)
            }
            res <- (mu * (1 - mu))/sigma^2 - 1
          }
          res
        })
    )

}


#### Niches ####

#' Imports data on water table depth growing niches for Sphagnum species
#'
#' @param file
#'
#' @export
leaching_data_get_niche_wtd <- function(file) {

  # water table depths where Sphagnum species grow
  res <-
    read.csv(file) %>%
    dplyr::select(-c(2, 4, 9)) %>%
    setNames(nm = c("taxon_rank_value", "pH", paste0(c("Ca", "Mg", "Na", "K"), "_relativa_mass"), "water_table_depth_1", "sample_size")) %>%
    dplyr::mutate(
      water_table_depth_1_error_type = "sd", #---note: assumed
      water_table_depth_1_error =
        water_table_depth_1 %>%
        stringr::str_extract("\\(\\d+\\.?\\d?\\)") %>%
        stringr::str_remove("\\(") %>%
        stringr::str_remove("\\)") %>%
        as.numeric(),
      water_table_depth_1 =
        water_table_depth_1 %>%
        stringr::str_remove(" \\(\\d+\\.?\\d?\\)") %>%
        as.numeric()
    ) %>%
    dplyr::rename(
      water_table_depth_to_surface_niche_johnson2015 = "water_table_depth_1",
      water_table_depth_to_surface_niche_johnson2015_error = "water_table_depth_1_error",
      water_table_depth_to_surface_niche_johnson2015_error_type = "water_table_depth_1_error_type"
    )

  res

}


#' Adds growing niche water table depths to the main table
#'
#' @param x Main table
#'
#' @param leaching_data_sphagnum_niches Output from `leaching_data_get_niche_wtd()`.
#'
#' @export
leaching_data_add_niche_wtd <- function(x, leaching_data_sphagnum_niches) {

  dplyr::left_join(
    x,
    leaching_data_sphagnum_niches %>%
      dplyr::select(taxon_rank_value, water_table_depth_to_surface_niche_johnson2015),
    by = "taxon_rank_value"
  )

}



#### Simulated data (for sensitivity analysis) ####

#' Simulates data for the sensitivity analysis
#'
#' @param id_dataset Character value specifying which dataset to simulate.
#'
#' @export
leaching_get_simulated_data <- function(id_dataset, mass_relative_mass_offset, x = NULL, x_model_draws = NULL) {

  res <-
  switch(
    id_dataset,
    "leaching_data_sim_1" = {

      expand.grid(
        k_2_sim = c(0.01, 0.05, 0.15),
        l_2_sim = c(0.01, 0.05, 0.15),
        alpha_2_sim = c(1.00001, 3),
        phi_2_sim = c(200),
        mass_relative_mass_precision_is_na = c(TRUE, FALSE),
        incubation_duration =
          list(
            `Design 1` = c(1, 2),
            `Design 2` = c(0.05, 1, 2),
            `Design 3` = c(0.05, 0.5, 1, 2, 3, 5)
          ),
        n_samples_per_sampling_date = c(5, 10)
      ) %>%
      dplyr::mutate(
        id_sample_incubation_start = seq_len(nrow(.)),
        incubation_duration_design = names(incubation_duration),
        taxon_rank_value =
          paste0(k_2_sim, "_", alpha_2_sim),
        id_dataset =
          paste0(phi_2_sim, "_", l_2_sim, "_", incubation_duration_design),
        id_citation = id_dataset,
        speciesxstudies = paste0(taxon_rank_value, "_", id_dataset)
      ) %>%
      tidyr::unnest(incubation_duration) %>%
      dplyr::mutate(
        mass_relative_mass_mu_sim =
          m70(
            t = incubation_duration,
            k_1 = k_2_sim,
            l_1 = l_2_sim,
            alpha_1 = alpha_2_sim,
            s = mass_relative_mass_offset
          ),
        mass_relative_mass =
          purrr::map(seq_len(nrow(.)), function(i) {
            leaching_simulate_beta_mu_phi_1(
              mu = mass_relative_mass_mu_sim[[i]],
              phi = phi_2_sim[[i]],
              n = n_samples_per_sampling_date[[i]]
            )
          }),
        incubation_duration = incubation_duration * 365
      ) %>%
      tidyr::unnest(mass_relative_mass) %>%
      dplyr::mutate(
        id_sample = seq_len(nrow(.)),
        mass_relative_mass_precision =
          dplyr::case_when(
            mass_relative_mass_precision_is_na ~ NA_real_,
            TRUE ~ mass_relative_mass_precision
          )
      )

    },
    "leaching_data_sim_2" = {

      index <- sample(nrow(posterior::draws_of(x_model_draws$k_2[[1]])), size = 1)

      x_model_draws %>%
        dplyr::select(-id_fit) %>%
        dplyr::filter(! is.na(m_rep)) %>%
        dplyr::select(id_sample, id_sample_incubation_start, taxon_rank_value, id_citation, speciesxstudies, k_2, l_2, alpha_2, phi_2_p2, mass_relative_mass_precision, incubation_duration, mass_relative_mass_sample_size) %>%
        dplyr::rename(
          k_2_sim = "k_2",
          phi_2_sim = "phi_2_p2",
          l_2_sim = "l_2",
          alpha_2_sim = "alpha_2"
        ) %>%
        dplyr::mutate(
          n_samples_per_sampling_date = ifelse(is.na(mass_relative_mass_sample_size), 5, mass_relative_mass_sample_size),
          mass_relative_mass_precision_is_na = ifelse(is.na(mass_relative_mass_precision), TRUE, FALSE),
          dplyr::across(
            dplyr::all_of(c("k_2_sim", "l_2_sim", "alpha_2_sim", "phi_2_sim")),
            function(.x) {
              posterior::draws_of(.x)[index, ]
            }
          ),
          mass_relative_mass_mu_sim =
            m70(
              t = incubation_duration/365,
              k_1 = k_2_sim,
              l_1 = l_2_sim,
              alpha_1 = alpha_2_sim,
              s = mass_relative_mass_offset
            ),
          mass_relative_mass =
            purrr::map(seq_len(nrow(.)), function(i) {
              leaching_simulate_beta_mu_phi_1(
                mu = mass_relative_mass_mu_sim[[i]],
                phi = phi_2_sim[[i]],
                n = n_samples_per_sampling_date[[i]]
              )
            })
        ) %>%
        dplyr::select(-mass_relative_mass_precision) %>%
        tidyr::unnest(mass_relative_mass) %>%
        dplyr::mutate(
          mass_relative_mass_precision =
            dplyr::case_when(
              mass_relative_mass_precision_is_na ~ NA_real_,
              TRUE ~ mass_relative_mass_precision
            )
        )

    }
  )

  res

}


#' Helper function to simulate from a beta distribution
#'
#' Simulates a sample average and precision parameter from a beta distrribution.
#'
#' @param mu
#'
#' @param phi
#'
#' @param n
#'
#' @export
leaching_simulate_beta_mu_phi_1 <- function(mu, phi, n) {

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  res <- rbeta(n, shape1, shape2)
  mu <- mean(res)
  sigma <- sd(res)
  if(sigma^2 >= (mu * (1 - mu))) {
    sigma <- sqrt(mu * (1 - mu)) * 0.99
  }
  precision <- (mu * (1 - mu))/sigma^2 - 1

  tibble::tibble(
    mass_relative_mass = mu,
    mass_relative_mass_precision = precision
  )

}



#### Combine ####

#' Combine previous functions
#'
#' @export
leaching_prepare_data <- function(file, mass_relative_mass_offset, leaching_data_sphagnum_niches) {

  res <-
    leaching_load_data_database(file) %>%
    leaching_filter_data_1() %>%
    leaching_add_metadata_1() %>%
    leaching_data_make_wide() %>%
    leaching_get_average_wtd() %>%
    leaching_filter_data_2() %>%
    leaching_data_adjust_taxon_rank_value() %>%
    leaching_data_add_niche_wtd(leaching_data_sphagnum_niches = leaching_data_sphagnum_niches) %>%
    leaching_prepare_remaining_masses(mass_relative_mass_offset)

  res

}
