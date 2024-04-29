#' Computes the water content of peat layers from their `water_table_depth_1` using the modified Granberg model [@Granberg.1999, @Kettridge.2007]
#'
#' @param layer_depth_midpoint_1 A numeric vector with the depth of the
#' midpoints of peat layers below the ground surface [cm].
#'
#' @param layer_total_porosity_1 A numeric vector with the total porosity of
#' peat layers [L L\eqn{^{-1}}].
#'
#' @param layer_water_table_depth_1 A vector with the distance between the
#' midpoint of a layer to the water table level [cm]. Negative values indicate
#' that the midpoint of the layer is below the water table level.
#'
#' @param minimum_water_content_at_surface_1 A numeric vector with the minimum
#' water content a peat layer at a depth of 0 cm can have [L L\eqn{^{-1}}].
#'
#' @return A numeric vector (`layer_water_content_1`) with the water content for
#' each peat layer [L L\eqn{^{-1}}].
#'
#' @examples
#' layer_depth_midpoint_1 = pdpm_sample_data_layers_from_1[, 1] + 0.5 * pdpm_sample_data_layers_from_1[, 3]
#' water_table_depth_to_surface_1 = 10
#' minimum_water_content_at_surface_1 = 0.2
#'
#' layer_water_content_1 <-
#'   m7(
#'     layer_depth_midpoint_1 = layer_depth_midpoint_1,
#'     layer_total_porosity_1 = pdpm_sample_data_layers_from_1[, 6],
#'     water_table_depth_to_surface_1 = water_table_depth_to_surface_1,
#'     minimum_water_content_at_surface_1 = minimum_water_content_at_surface_1
#'   )
#'
#' @export
m7 <- function(layer_depth_midpoint_1, layer_total_porosity_1, water_table_depth_to_surface_1, minimum_water_content_at_surface_1) {

  J = length(layer_depth_midpoint_1)
  water_content_at_surface_1 = numeric(1L)
  res = numeric(J)
  water_table_depth_to_surface_1 = water_table_depth_to_surface_1/100
  layer_depth_midpoint_1 = layer_depth_midpoint_1/100

  # compute water content at surface
  if(water_table_depth_to_surface_1 < 0) {
    water_content_at_surface_1 = 1.0
  } else {
    water_content_at_surface_1 = 0.15 * water_table_depth_to_surface_1^(-0.28) #---note: assumes that the first wtd value corresponds to the wtd of the surface layer
  }
  if(water_content_at_surface_1 < minimum_water_content_at_surface_1) {
    water_content_at_surface_1 = minimum_water_content_at_surface_1
  }
  if(water_content_at_surface_1 > layer_total_porosity_1[[1]]) {
    water_content_at_surface_1 = layer_total_porosity_1[[1]]
  }

  # compute water content
  res = water_content_at_surface_1 + (layer_total_porosity_1 - water_content_at_surface_1) * (layer_depth_midpoint_1 / water_table_depth_to_surface_1)^2
  for(j in 1:J) {
    if(res[j] > layer_total_porosity_1[[j]] || water_table_depth_to_surface_1 <= 0) {
      res[[j]] = layer_total_porosity_1[j]
    }
  }

  res
}


#' Same as `m7`, but optimized for rvars
#'
#' Also allows to set a separate surface peat porosity for each sample.
#'
#' @export
m7_rvars <- function(layer_depth_midpoint_1, layer_total_porosity_1, water_table_depth_to_surface_1, minimum_water_content_at_surface_1, layer_total_porosity_at_surface) {

  water_table_depth_to_surface_1 = water_table_depth_to_surface_1/100
  layer_depth_midpoint_1 = layer_depth_midpoint_1/100

  # compute water content at surface
  water_content_at_surface_1 <-
    posterior::rvar_ifelse(
      water_table_depth_to_surface_1 < 0,
      1.0,
      0.15 * water_table_depth_to_surface_1^(-0.28) #---note: assumes that the first wtd value corresponds to the wtd of the surface layer
    )

  water_content_at_surface_1 <-
    posterior::rvar_ifelse(
      water_content_at_surface_1 < minimum_water_content_at_surface_1,
      minimum_water_content_at_surface_1,
      water_content_at_surface_1
    )

  water_content_at_surface_1 <-
    posterior::rvar_ifelse(
      water_content_at_surface_1 > layer_total_porosity_at_surface,
      layer_total_porosity_at_surface,
      water_content_at_surface_1
    )

  # compute water content
  res = water_content_at_surface_1 + (layer_total_porosity_1 - water_content_at_surface_1) * (layer_depth_midpoint_1 / water_table_depth_to_surface_1)^2
  res <-
    posterior::rvar_ifelse(
      res > layer_total_porosity_1 | water_table_depth_to_surface_1 <= 0,
      layer_total_porosity_1,
      res
    )

  res
}



#' The decomposition model in the Holocene Peatland Model [@Frolking.2010].
#'
#' @inheritParams m69
#'
#' @param layer_water_table_depth_1
#'
#' @param m68_p1
#'
#' @param m68_p2
#'
#' @param m68_p3
#'
#' @param m68_p4
#'
#' @param m68_p5
#'
#' @return `layer_decomposition_rate_1`.
#'
#' @examples
#' layer_decomposition_rate_1 <-
#'   m68(
#'     layer_degree_of_saturation_1 = c(0, 0.2, 0.45, 0.7, 1, 1, 1, 1),
#'     layer_water_table_depth_1 = c(80, 40, 20, 10, 0, -10, -50, -100),
#'     m69_p1 = 0.45,
#'     m69_p2 = 2.31,
#'     m68_p1 = 0.001,
#'     m68_p2 = 0.3,
#'     m68_p3 = 0.08,
#'     m68_p4 = 1,
#'     m68_p5 = 1
#'   )
#'
#' @export
m68 <- function(layer_degree_of_saturation_1, layer_water_table_depth_1, m69_p1 = 0.45, m69_p2 = 2.31, m68_p1 = 0.001, m68_p2 = 0.3, m68_p3, m68_p4, m68_p5) {

  res <-
    m69(
      layer_degree_of_saturation_1 = layer_degree_of_saturation_1,
      m69_p1 = m69_p1,
      m69_p2 = m69_p2
    )

  res1 <-
    m69(
      layer_degree_of_saturation_1 = 1.0,
      m69_p1 = m69_p1,
      m69_p2 = m69_p2
    )

  res <-
    posterior::rvar_ifelse(
      layer_water_table_depth_1 <= 0,
      (m68_p1 + (res1 - m68_p1) * exp(layer_water_table_depth_1/100/m68_p2)),
      res
    )

  res <- m68_p3 * (m68_p4/m68_p5) * res

  res

}



#' The part of the decomposition model in the Holocene Peatland Model [@Frolking.2010] active above the water table.
#'
#' @param layer_degree_of_saturation_1
#'
#' @param m69_p1
#'
#' @param m69_p2
#'
#' @return `layer_decomposition_rate_1`.
#'
#' @examples
#' layer_decomposition_rate_1 <-
#'   m69(
#'     layer_degree_of_saturation_1 = c(0, 0.2, 0.45, 0.7, 1, 1, 1, 1),
#'     m69_p1 = 0.45,
#'     m69_p2 = 2.31
#'   )
#'
#' @export
m69 <- function(layer_degree_of_saturation_1, m69_p1 = 0.45, m69_p2 = 2.31) {

  1 - m69_p2 * (layer_degree_of_saturation_1 - m69_p1)^2

}



