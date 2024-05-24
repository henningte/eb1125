vector m70(vector t, vector k_1, vector l_1, vector alpha_1, real s) {

  int N = size(t);
  vector[N] res;

  res = (1.0 - l_1) ./ (1.0 + (alpha_1 - 1.0) .* k_1 .* t)^(1.0 ./ (alpha_1 - 1.0));
  for(n in 1:N) {
    if(res[n] >= 1.0) {
      res[n] = 1.0 - s;
    }
  }

  return(res);

}
