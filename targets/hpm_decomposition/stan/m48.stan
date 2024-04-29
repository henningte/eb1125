// Computes the water content of peat layers from their `water_table_depth_1` using the modified Granberg model [@Granberg.1999, @Kettridge.2007]
// * `layer_depth_midpoint_1`: A vector with an element for each peat layer for which to
//   compute the moisture content. Elements are the depth (for a peat
//   layer, one can use the average peat depth) [cm].
// * `layer_total_porosity_1`: A vector with an element for each peat layer for which to
//   compute the moisture content. Elements are the total porosity [L/L].
// * `layer_water_table_depth_1`: A vector with an element for each peat layer for which to
//   compute the moisture content. Elements are the water table depth of the
//   respective peat layer [cm].
// * `minimum_water_content_at_surface_1`: A numeric value representing the minimum
//   water content [L/L] at the peat surface.
vector m48(vector layer_depth_midpoint_1, vector layer_water_table_depth_1, vector layer_total_porosity_1, vector layer_total_porosity_1_at_surface, vector minimum_water_content_at_surface_1) {

  int J = size(layer_depth_midpoint_1);
  vector[J] water_content_at_surface_1;
  vector[J] res;
  vector[J] wtd2 = layer_water_table_depth_1/100;

  // compute water_content_at_surface_1
  for(j in 1:J) {

    if(wtd2[j] <= 0) { //---note: water table is above peat surface
      res[j] = layer_total_porosity_1[j];
    } else {

      water_content_at_surface_1[j] = 0.15 * wtd2[j]^(-0.28);

      if(water_content_at_surface_1[j] < minimum_water_content_at_surface_1[j]) {
        water_content_at_surface_1[j] = minimum_water_content_at_surface_1[j];
      }
      if(water_content_at_surface_1[j] > layer_total_porosity_1_at_surface[j]) {
        water_content_at_surface_1[j] = layer_total_porosity_1_at_surface[j];
      }

      res[j] = water_content_at_surface_1[j] + (layer_total_porosity_1_at_surface[j] - water_content_at_surface_1[j]) * ((layer_depth_midpoint_1[j]/100)./wtd2[j])^2;
      if(res[j] > layer_total_porosity_1[j]) {
        res[j] = layer_total_porosity_1[j];
      }

    }

  }

  return(res);
}
