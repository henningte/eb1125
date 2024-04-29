//// The decomposition model in the Holocene Peatland Model [@Frolking.2010]. Takes as input `layer_degree_of_saturation_1`, `layer_water_table_depth_1`, and parameters.
//
// @inheritsParams m69
vector m68(
  vector layer_degree_of_saturation_1,
  vector layer_water_table_depth_1,
  real m69_p1,
  real m69_p2,
  real m68_p1,
  real m68_p2,
  vector m68_p3
) {

 int N = size(layer_degree_of_saturation_1);
 vector[N] res;

 res =
   m69(
     layer_degree_of_saturation_1,
     m69_p1,
     m69_p2
   );

  for(n in 1:N) {
    if(layer_water_table_depth_1[n] <= 0) {
      res[n] = m68_p1 + (res[n] - m68_p1) * exp(layer_water_table_depth_1[n]/100/m68_p2);
    }
  }

  res = m68_p3 .* res;

  return(res);

}

