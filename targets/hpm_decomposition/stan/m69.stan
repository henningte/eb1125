//// The part of the decomposition model in the Holocene Peatland Model [@Frolking.2010] active above the water table.
//
vector m69(
  vector layer_degree_of_saturation_1,
  real m69_p1,
  real m69_p2
) {

 int N = size(layer_degree_of_saturation_1);
 vector[N] res = 1.0 - m69_p2 * (layer_degree_of_saturation_1 - m69_p1)^2;

 return(res);

}

