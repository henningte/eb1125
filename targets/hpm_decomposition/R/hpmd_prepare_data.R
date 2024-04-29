#### HPM microhabitats ####

#' Adds growing niche water table depth errors to the main table
#'
#' @param x Main table
#'
#' @param leaching_data_sphagnum_niches
#'
#' @export
hpmd_data_add_niche_wtd <- function(x, leaching_data_sphagnum_niches) {

  dplyr::left_join(
    x,
    leaching_data_sphagnum_niches %>%
      dplyr::select(taxon_rank_value, water_table_depth_to_surface_niche_johnson2015_error, water_table_depth_to_surface_niche_johnson2015_error_type),
    by = "taxon_rank_value"
  )

}


#' Defines Holocene Peatland Model microhabitats
#'
#' Based on niche water table depths.
#'
#' @param leaching_data_sphagnum_niches `leaching_data_sphagnum_niches`.
#'
#' @return `leaching_data_sphagnum_niches` with selected columns and an
#' assignment of species to microhabitats.
#'
#' @export
hpmd_define_hpm_microhabitats <- function(leaching_data_sphagnum_niches = leaching_data_sphagnum_niches) {

  res <-
    leaching_data_sphagnum_niches %>%
    dplyr::select(taxon_rank_value, water_table_depth_to_surface_niche_johnson2015, water_table_depth_to_surface_niche_johnson2015_error) %>%
    dplyr::mutate(
      smaller5 = pnorm(q = 5, mean = water_table_depth_to_surface_niche_johnson2015, sd = water_table_depth_to_surface_niche_johnson2015_error) * 100,
      smaller15larger5 = pnorm(q = 15, mean = water_table_depth_to_surface_niche_johnson2015, sd = water_table_depth_to_surface_niche_johnson2015_error) * 100 - smaller5,
      larger15 = 100 - pnorm(q = 15, mean = water_table_depth_to_surface_niche_johnson2015, sd = water_table_depth_to_surface_niche_johnson2015_error) * 100,
      dplyr::across(
        dplyr::starts_with(c("smaller", "larger")),
        \(.x) round(.x, digits = 1)
      ),
      hpm_microhabitat =
        dplyr::case_when(
          larger15 > 15 & smaller5 > 15 & smaller15larger5 > 15 ~ "Lawn",
          smaller5 > smaller15larger5 & smaller5 > larger15 ~ "Hollow",
          smaller15larger5 > smaller5 & smaller15larger5 > larger15 ~ "Lawn",
          larger15 > smaller5 & larger15 > smaller15larger5 ~ "Hummock"
        )
    )

  res

}


#' Adds hpm_microhabitat to the main table
#'
#' @param x
#'
#' @param leaching_data_hpm_microhabitat
#'
#' @export
hpmd_data_add_hpm_microhabitat <- function(x, hpmd_data_hpm_microhabitat) {

  # assign samples to standard HPM microhabitats
  res <-
    x %>%
    dplyr::left_join(
      hpmd_data_hpm_microhabitat %>%
        dplyr::select(taxon_rank_value, hpm_microhabitat),
      by = "taxon_rank_value"
    ) %>%
    dplyr::mutate(
      hpm_microhabitat =
        dplyr::case_when(
          taxon_rank_value == "Sphagnum" ~ "Hummock",
          TRUE ~ hpm_microhabitat
        ) %>%
        factor(levels = c("Hollow", "Lawn", "Hummock"))
    )

  # assign samples to custom hpm microhabitats (each species is an own HPM micrhabitat class)
  res <-
    res %>%
    dplyr::mutate(
      hpm_microhabitat2 =
        dplyr::case_when(
          taxon_rank_value == "Sphagnum" & sample_depth_lower == 10 ~ "Sphagnum_1",
          taxon_rank_value == "Sphagnum" & sample_depth_lower == 30 ~ "Sphagnum_2",
          TRUE ~ taxon_rank_value
        ),
      hpm_microhabitat2 =
        dplyr::case_when(
          index_hpm ~ hpm_microhabitat2,
          TRUE ~ NA_character_
        )
    )

  res

}


#### Indices ####

#' Identifies samples which can be compared to the HPM
#'
#' @param x The main data table
#'
#' @export
hpmd_identify_samples_for_hpm <- function(x, leaching_data_hpm_microhabitat) {

  x %>%
    dplyr::mutate(
      index_hpm = ! is.na(water_table_depth_to_surface_average)
    )

}



#' Identify samples eligible for cross-validation
#'
#' @param x The main data table.
#'
#' @export
hpmd_add_cross_validation_indices <- function(x) {

  x %>%
    dplyr::left_join(
      x %>%
        dplyr::filter(index_hpm) %>%
        dplyr::group_by(hpm_microhabitat2) %>%
        dplyr::summarise(
          hpm_n_study_for_hpm_microhabitat2 = length(unique(id_citation)),
          .groups = "drop"
        ),
      by = "hpm_microhabitat2"
    ) %>%
    dplyr::mutate(
      index_hpm_cross_validation = index_hpm & hpm_n_study_for_hpm_microhabitat2 > 1, # measurements which are principally eligible for cross-validation
      hpm_cross_validation_id_block =
        id_citation %>%
        factor(levels = unique(id_citation[index_hpm_cross_validation])),
      hpm_cross_validation_id_block =
        ifelse(! index_hpm_cross_validation, NA, hpm_cross_validation_id_block)
    ) %>%
    dplyr::select(-hpm_n_study_for_hpm_microhabitat2)

}




#### Combine ####


#' Prepares the main data table
#'
#' @param x Corresponds to target `leaching_data_from_database` from folder
#' `leaching`.
#'
#' @export
hpmd_prepare_data <- function(x, hpmd_data_hpm_microhabitat = hpmd_data_hpm_microhabitat) {

  x %>%
    hpmd_identify_samples_for_hpm() %>%
    hpmd_data_add_hpm_microhabitat(hpmd_data_hpm_microhabitat = hpmd_data_hpm_microhabitat) %>%
    hpmd_add_cross_validation_indices()

}

