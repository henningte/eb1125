#' Returns a data frame with standard HPM parameter values
#'
#' @export
hpmd_get_hpm_standard_parameter_values <- function() {

  tibble::tibble(
    variable = c("m69_p1","m69_p2", "m68_p1", "m68_p2", rep("m68_p3", 3)),
    hpm_microhabitat = c(rep(NA_character_, 4), c("Hollow", "Lawn", "Hummock"))
  ) %>%
    dplyr::mutate(
      value =
        dplyr::case_when(
          variable == "m69_p1" ~ 0.45,
          variable == "m69_p2" ~ 2.31,
          variable == "m68_p1" ~ 0.001,
          variable == "m68_p2" ~ 0.3,
          variable == "m68_p3" & hpm_microhabitat == "Hollow" ~ 0.13,
          variable == "m68_p3" & hpm_microhabitat == "Lawn" ~ 0.08,
          variable == "m68_p3" & hpm_microhabitat == "Hummock" ~ 0.06
        ),
      variable_pretty =
        dplyr::case_when(
          variable == "m69_p1" ~ "<i>W<sub>opt</sub></i> (L<sub>water</sub> L<sup>-1</sup><sub>pores</sub>)",
          variable == "m69_p2" ~ "<i>c<sub>1</sub></i> (-)",
          variable == "m68_p1" ~ "<i>f<sub>min</sub></i> (yr<sup>-1</sup>)",
          variable == "m68_p2" ~ "<i>c<sub>2</sub></i> (m)",
          variable == "m68_p3" & hpm_microhabitat == "Hollow" ~ "<i>k<sub>0, hollow</sub></i> (yr<sup>-1</sup>)",
          variable == "m68_p3" & hpm_microhabitat == "Lawn" ~ "<i>k<sub>0, lawn</sub></i> (yr<sup>-1</sup>)",
          variable == "m68_p3" & hpm_microhabitat == "Hummock" ~ "<i>k<sub>0, hummock</sub></i> (yr<sup>-1</sup>)"
        )
    )

}


#' Returns a data frame with peat properties assumed here
#'
#' @export
hpmd_get_peat_properties <- function() {

  tibble::tibble(
    layer_total_porosity_1_average = 0.8,
    layer_total_porosity_1_error = 0.1,
    layer_total_porosity_1_error_type = "sd",
    minimum_water_content_at_surface_1_average = 0.05,
    minimum_water_content_at_surface_1_error = 0.05, #---todo
    minimum_water_content_at_surface_1_error_type = "sd",
    water_table_depth_to_surface_1_error = 3,
    water_table_depth_to_surface_1_error_type = "sd"
  )

}


#' Converts model ids into model names
#'
#' @export
hpmd_model_id_to_name <- function(x) {

  id_models <- paste0(3, "-", 1:5)

  purrr::map_chr(x, function(.x) {
    index <- id_models == .x
    if(sum(index) == 0) {
      NA_character_
    } else {
      c("HPMf", "HPMf-LE-peat", "HPMe-LE-peat", "HPMe-LE-peat-l0", "HPMe-LE-peat-l0-outlier")[index]
    }

  })

}
